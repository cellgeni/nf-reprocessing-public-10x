package cellgeni.lineagekeep

import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.CountDownLatch
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.TimeUnit
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session
import nextflow.dag.DAG
import nextflow.extension.GroupKey
import nextflow.processor.ErrorStrategy
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ProcessConfigV1
import nextflow.script.ProcessConfigV2
import nextflow.script.params.FileOutParam
import nextflow.script.params.OutParam
import nextflow.trace.TraceObserverV2
import nextflow.trace.event.FilePublishEvent
import nextflow.trace.event.TaskEvent
import nextflow.util.Threads

/**
 * Deletes the work directories of tasks whose downstream lineage succeeded, keeps every
 * task upstream of a final failure, and at the end writes the failures and a re-run table.
 *
 * Observer callbacks arrive on several Nextflow threads. Each one only snapshots what it
 * needs from the task and queues it; a single worker thread owns the {@link TaskGraph} and
 * does the deleting, so neither needs locking and no Nextflow thread waits on Lustre.
 *
 * The DAG walk in {@link #processOutputs} follows nf-boost's CleanupObserver
 * (Ben Sherman, Apache-2.0): https://github.com/bentsherman/nf-boost
 */
@Slf4j
@CompileStatic
class LineageKeepObserver implements TraceObserverV2 {

    private final LineageKeepConfig config
    private Session session
    private String workRoot
    private Path outDir
    private TaskGraph graph
    private WorkDirCleaner cleaner
    private final Map<String,ProcessOutputs> processes = new HashMap<>()
    private final List<Tuple2<Pattern,String>> scratch = []

    private final LinkedBlockingQueue<Runnable> events = new LinkedBlockingQueue<>()
    private final CountDownLatch done = new CountDownLatch(1)
    private Thread worker
    private volatile boolean broken

    LineageKeepObserver(LineageKeepConfig config) {
        this.config = config
    }

    /** Per output index of a process: who consumes it, and whether it may be published. */
    @PackageScope
    static class ProcessOutputs {
        final List<Set<String>> consumers
        final List<Boolean> publishable
        ProcessOutputs(List<Set<String>> consumers, List<Boolean> publishable) {
            this.consumers = consumers
            this.publishable = publishable
        }
    }

    /** What the worker needs from a task, read on the calling thread. */
    @PackageScope
    static class TaskInfo {
        String key
        String process
        String name
        String hash
        String workDir
        Set<String> inputDirs = new LinkedHashSet<>()
        /** output path -> processes that read it */
        Map<String,Set<String>> outputs = new LinkedHashMap<>()
        List<String> publishable = []
        String dataset
        String sample
        List<String> sampleList = []
        boolean success
        boolean retry
        String exit
    }

    // ---- lifecycle -------------------------------------------------------

    @Override
    void onFlowCreate(Session session) {
        this.session = session
        this.workRoot = normalize(session.workDir)
        this.outDir = config.outputDir
            ? Path.of(config.outputDir).toAbsolutePath()
            : session.outputDir.resolve('cleanup').toAbsolutePath()
        this.graph = new TaskGraph(config.barrierProcesses, config.singleUseProcesses)
        this.cleaner = new WorkDirCleaner(outDir.resolve('deleted.tsv'), config.dryRun, config.keepLogs)
        for( String spec : config.scratch ) {
            final i = spec.indexOf(':')
            if( i <= 0 ) {
                log.warn "[nf-lineage-keep] ignoring scratch entry `${spec}` -- expected PROCESS:glob"
                continue
            }
            scratch << new Tuple2<Pattern,String>(LineageKeepConfig.wildcard(spec.substring(0, i)), spec.substring(i + 1))
        }
        log.info "[nf-lineage-keep] ${config.dryRun ? 'DRY RUN -- nothing will be deleted' : 'deleting'}; log: ${outDir.resolve('deleted.tsv')}"
    }

    @Override
    void onFlowBegin() {
        try {
            analyseDag(session.dag, session.outputs.values())
            writeProcessTable()
        }
        catch( Throwable e ) {
            // without the consumer map nothing can be judged safe, so the plugin stands down
            broken = true
            log.warn "[nf-lineage-keep] unable to read the process graph, no work directory will be deleted -- ${e.message ?: e}", e
            return
        }
        worker = Threads.start('nf-lineage-keep') {
            try {
                while( !done.await(config.interval.toMillis(), TimeUnit.MILLISECONDS) )
                    drainAndClean(false)
            }
            catch( Throwable e ) {
                broken = true
                log.warn "[nf-lineage-keep] stopped after an unexpected error, no further work directory will be deleted -- ${e.message ?: e}", e
            }
        }
    }

    @Override
    void onFlowComplete() {
        if( worker == null )
            return
        done.countDown()
        worker.join()
        try {
            // every process has ended, so publishing is over too; but an aborted run stopped
            // its processes half-way and nothing about it can be called finished
            final complete = session.isSuccess() && !session.isAborted() && !session.isCancelled()
            drainAndClean(complete)
            writeReports(complete)
        }
        catch( Throwable e ) {
            log.warn "[nf-lineage-keep] final pass failed -- ${e.message ?: e}", e
        }
        finally {
            cleaner.close()
        }
    }

    // ---- events (Nextflow threads) ---------------------------------------

    @Override
    void onTaskPending(TaskEvent event) {
        if( broken ) return
        final info = snapshot(event.handler.task, false)
        events.put({ register(info) } as Runnable)
    }

    @Override
    void onTaskCached(TaskEvent event) {
        if( broken ) return
        final info = snapshot(event.handler.task, true)
        info.success = true
        events.put({ register(info); complete(info) } as Runnable)
    }

    @Override
    void onTaskComplete(TaskEvent event) {
        if( broken ) return
        final task = event.handler.task
        final info = snapshot(task, true)
        info.success = !task.failed && !task.aborted && task.isSuccess()
        info.retry = task.errorAction == ErrorStrategy.RETRY
        info.exit = task.aborted ? 'aborted' : String.valueOf(task.exitStatus)
        events.put({ complete(info) } as Runnable)
    }

    @Override
    void onFilePublish(FilePublishEvent event) {
        if( broken || event.source == null ) return
        final path = normalize(event.source)
        events.put({ graph.published(path) } as Runnable)
    }

    @Override
    void onProcessTerminate(TaskProcessor process) {
        if( broken ) return
        final name = process.name
        events.put({ graph.terminated(name) } as Runnable)
    }

    // ---- worker thread ---------------------------------------------------

    private void register(TaskInfo info) {
        final n = graph.register(info.key, info.process, info.workDir, info.inputDirs,
                                 info.name, info.hash, info.dataset, info.sample)
        if( info.sampleList )
            n.sampleList = info.sampleList
    }

    private void complete(TaskInfo info) {
        purgeScratch(info)
        if( info.success )
            graph.succeeded(info.key, info.outputs, info.publishable)
        else if( !info.retry )
            graph.failed(info.key, info.exit)
    }

    private void purgeScratch(TaskInfo info) {
        if( !scratch || info.workDir == null ) return
        final name = TaskGraph.simpleName(info.process)
        final globs = scratch.findAll { Tuple2<Pattern,String> s -> s.v1.matcher(name).matches() }.collect { it.v2 }
        if( globs )
            cleaner.purge(Path.of(info.workDir), globs, info.process, info.name, info.hash)
    }

    private void drainAndClean(boolean flowComplete) {
        final List<Runnable> batch = []
        events.drainTo(batch)
        for( Runnable r : batch )
            r.run()
        if( broken )
            return
        for( TaskGraph.Cleanup c : graph.collectDeletable(flowComplete) ) {
            final n = c.node
            if( !c.bulk ) {
                cleaner.remove(c.paths, n.process, n.name, n.hash)
                continue
            }
            final dirs = config.keepFailedAttempts ? [n.workDir] : n.attemptDirs
            for( String d : dirs )
                cleaner.clean(Path.of(d), c.paths, n.process, n.name, n.hash)
        }
        for( TaskGraph.Node n : graph.pinnedAfterDelete )
            log.warn "[nf-lineage-keep] ${n.name} was deleted before a failure downstream pinned it -- check singleUseProcesses for ${n.group}"
        graph.pinnedAfterDelete.clear()
    }

    // ---- task snapshot ---------------------------------------------------

    private TaskInfo snapshot(TaskRun task, boolean withOutputs) {
        final info = new TaskInfo()
        info.process = task.processor.name
        info.key = "${info.process}/${task.index}".toString()
        info.name = task.name
        info.hash = task.hash?.toString()
        info.workDir = task.workDir ? normalize(task.workDir) : null
        for( def holder : task.inputFiles ) {
            addInputDir(info.inputDirs, holder.storePath)
            addInputDir(info.inputDirs, holder.sourcePath)
        }
        readKeys(task, info)
        if( withOutputs )
            readOutputs(task, info)
        return info
    }

    private void addInputDir(Set<String> dirs, Path path) {
        if( path == null ) return
        final dir = taskDirOf(normalize(path))
        if( dir != null ) dirs << dir
    }

    /**
     * `<workRoot>/ab/cdef…/anything` -> `<workRoot>/ab/cdef…`; null outside the work dir.
     */
    @PackageScope
    static String taskDirOf(String path, String workRoot) {
        if( !path.startsWith(workRoot + '/') ) return null
        final parts = path.substring(workRoot.length() + 1).split('/', 3)
        return parts.length >= 2 ? "${workRoot}/${parts[0]}/${parts[1]}".toString() : null
    }

    private String taskDirOf(String path) { taskDirOf(path, workRoot) }

    private void readOutputs(TaskRun task, TaskInfo info) {
        final state = processes[info.process]
        for( Map.Entry<FileOutParam,Object> e : task.getOutputsByType(FileOutParam).entrySet() ) {
            final files = asPaths(e.value)
            if( !files ) continue
            final idx = e.key.index
            final known = state != null && idx < state.consumers.size()
            // an output the DAG walk did not see: assume an unknown reader, so it is never freed
            final readers = known ? state.consumers[idx] : (['?unknown'] as Set<String>)
            for( Path f : files ) {
                final path = normalize(f)
                info.outputs.computeIfAbsent(path) { new LinkedHashSet<String>() }.addAll(readers)
                if( known && state.publishable[idx] )
                    info.publishable << path
            }
        }
    }

    private static List<Path> asPaths(Object value) {
        if( value instanceof Path ) return [(Path)value]
        if( value instanceof Collection ) return ((Collection)value).findAll { it instanceof Path } as List<Path>
        return []
    }

    private void readKeys(TaskRun task, TaskInfo info) {
        try {
            final ctx = task.context
            if( ctx == null ) return
            final meta = unwrap(ctx.containsKey(config.metaVariable) ? ctx.get(config.metaVariable) : null)
            if( meta instanceof Map ) {
                final m = (Map)meta
                final dataset = str(m.get(config.datasetKey))
                final sample = str(m.get(config.sampleKey))
                final id = str(m.get(config.idKey))
                info.dataset = dataset ?: id
                info.sample = sample ?: (dataset ? id : null)
            }
            final list = unwrap(ctx.containsKey(config.sampleListVariable) ? ctx.get(config.sampleListVariable) : null)
            if( list instanceof Collection )
                info.sampleList = ((Collection)list).collect { it.toString().trim() }.findAll { it } as List<String>
            else if( list != null )
                info.sampleList = list.toString().tokenize(',')*.trim().findAll { it } as List<String>
        }
        catch( Exception e ) {
            log.debug "[nf-lineage-keep] no re-run key for ${task.name} -- ${e.message}"
        }
    }

    private static Object unwrap(Object o) {
        while( o instanceof GroupKey ) o = ((GroupKey)o).getGroupTarget()
        return o
    }

    private static String str(Object o) {
        o = unwrap(o)
        return o == null ? null : o.toString()
    }

    private static String normalize(Path p) {
        p.toAbsolutePath().normalize().toString()
    }

    // ---- DAG ---------------------------------------------------------------

    /**
     * For each process output, the processes that read it (through any chain of operators)
     * and whether it can reach a published channel.
     */
    private void analyseDag(DAG dag, Collection publishedChannels) {
        final Map<DAG.Vertex,List<DAG.Vertex>> successors = [:]
        final Map<DataflowWriteChannel,List<DAG.Edge>> edgesFrom = [:]
        for( DAG.Edge edge : dag.edges ) {
            successors.computeIfAbsent(edge.from) { [] }.add(edge.to)
            if( edge.channel instanceof DataflowWriteChannel )
                edgesFrom.computeIfAbsent((DataflowWriteChannel)edge.channel) { [] }.add(edge)
        }
        final published = new HashSet(publishedChannels)

        for( DAG.Vertex vertex : dag.vertices ) {
            final process = vertex.process
            if( process == null ) continue
            final outs = processOutputs(process)
            final List<Set<String>> consumers = []
            final List<Boolean> publishable = []
            for( OutParam param : outs ) {
                final ch = param.getOutChannel()
                // processes reached from this channel through operators
                final Set<String> found = new LinkedHashSet<>()
                final Set<DAG.Vertex> seen = new HashSet<>()
                final queue = new ArrayDeque<DAG.Vertex>()
                for( DAG.Edge e : edgesFrom.getOrDefault(ch, []) ) if( e.to != null ) queue << e.to
                while( queue ) {
                    final w = queue.poll()
                    if( !seen.add(w) ) continue
                    if( w.process != null ) found << w.process.name
                    else successors.getOrDefault(w, []).each { DAG.Vertex x -> if( x != null ) queue << x }
                }
                consumers << found
                publishable << reachesPublished(ch, edgesFrom, published)
            }
            processes[process.name] = new ProcessOutputs(consumers, publishable)
        }
    }

    private static boolean reachesPublished(DataflowWriteChannel start, Map<DataflowWriteChannel,List<DAG.Edge>> edgesFrom, Set published) {
        final seen = new HashSet()
        final queue = new ArrayDeque<DataflowWriteChannel>()
        queue << start
        while( queue ) {
            final ch = queue.poll()
            if( !seen.add(ch) ) continue
            if( ch in published ) return true
            for( DAG.Edge e : edgesFrom.getOrDefault(ch, []) ) {
                final w = e.to
                if( w == null || w.process != null ) continue
                for( def op : w.operators )
                    for( def out : op.getOutputs() )
                        if( out instanceof DataflowWriteChannel ) queue << (DataflowWriteChannel)out
            }
        }
        return false
    }

    private static List<OutParam> processOutputs(TaskProcessor process) {
        final config = process.config
        if( config instanceof ProcessConfigV2 )
            return ((ProcessConfigV2)config).getOutputs().getParams() as List<OutParam>
        if( config instanceof ProcessConfigV1 )
            return ((ProcessConfigV1)config).getOutputs() as List<OutParam>
        return []
    }

    // ---- reports -----------------------------------------------------------

    private void writeProcessTable() {
        Files.createDirectories(outDir)
        final sb = new StringBuilder('process\tgroup\toutput\tconsumers\tpublishable\n')
        processes.keySet().sort().each { String p ->
            final s = processes[p]
            for( int i = 0; i < s.consumers.size(); i++ )
                sb << p << '\t' << graph.groupOf(p) << '\t' << i << '\t' << s.consumers[i].join(',') << '\t' << s.publishable[i] << '\n'
        }
        Files.writeString(outDir.resolve('processes.tsv'), sb.toString())
    }

    private void writeReports(boolean complete) {
        Files.createDirectories(outDir)
        final failures = new StringBuilder('process\tname\thash\texit\tdataset\tsample\tworkdir\tkept_upstream\n')
        for( TaskGraph.Node n : graph.finalFailures() ) {
            final up = TaskGraph.ancestors(n).collect { TaskGraph.Node a -> a.workDir }.findAll { it }
            failures << [n.process, n.name, n.hash, n.exit, n.dataset, n.sample, n.workDir, up.join(',')]
                .collect { (it ?: '').toString().replace('\t', ' ') }.join('\t') << '\n'
        }
        Files.writeString(outDir.resolve('failed_tasks.tsv'), failures.toString())

        if( !complete ) {
            log.warn "[nf-lineage-keep] the run did not complete: rerun.tsv not written, as unfinished items would be missing from it"
            return
        }
        // the pipeline's own input format: tab-separated, one row per dataset,
        // sample_id comma-joined, and empty for the whole dataset
        final rerun = new StringBuilder('dataset_id\tsample_id\n')
        graph.rerunItems().each { String dataset, Set<String> samples ->
            rerun << dataset << '\t' << ('' in samples ? '' : samples.join(',')) << '\n'
        }
        Files.writeString(outDir.resolve('rerun.tsv'), rerun.toString())

        final pinned = graph.nodes.count { TaskGraph.Node n -> n.pinned }
        final GB = 1024L * 1024 * 1024
        log.info "[nf-lineage-keep] ${config.dryRun ? 'would have cleaned' : 'cleaned'} ${cleaner.dirs} work dirs, " +
            "${cleaner.totalBytes.intdiv(GB)} GB apparent / ${cleaner.totalFreed.intdiv(GB)} GB freed; " +
            "${graph.finalFailures().size()} failed tasks, ${pinned} tasks kept for them; re-run table: ${outDir.resolve('rerun.tsv')}"
    }
}
