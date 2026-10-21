package cellgeni.lineagekeep

import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j

/**
 * The run's task graph, as far as deletion is concerned. Pure bookkeeping, with no
 * Nextflow types in it, so that it can be unit-tested without a session.
 *
 * A node is one logical task: all attempts of `process/index`. Edges come from files: a task
 * whose input lives in another task's work directory consumes that task. Process-level
 * consumers (which processes read a process's output channels, from the DAG) are supplied
 * per node by the observer, because a consumer task may not exist yet when its producer
 * finishes.
 *
 * A node is OK when it succeeded, everything it emits that may be published has been
 * published, and every consumer it has or may still get, outside the barrier processes,
 * is itself OK. An OK node's work directory is cleaned, except for the declared outputs that
 * a consumer may still read (in practice: what the end-of-run barrier processes collect,
 * such as versions.yml or a STARsolo directory). Those are held, and deleted one by one
 * once nothing can read them any more. A final failure pins the failed task, its ancestors,
 * and every task that shares its re-run key; pinned nodes are never touched.
 */
@Slf4j
@CompileStatic
class TaskGraph {

    enum State { PENDING, SUCCEEDED, FAILED }

    static class Node {
        final String key
        final String process
        final String group
        State state = State.PENDING
        String name
        String hash
        String exit
        /** Work directories of every attempt, oldest first; the last one is current. */
        final List<String> attemptDirs = []
        final Set<Node> producers = new LinkedHashSet<>()
        final Set<Node> consumers = new LinkedHashSet<>()
        /** Output path -> processes that read it, for the outputs this task produced. */
        Map<String,Set<String>> outputs = Collections.<String,Set<String>>emptyMap()
        /** Processes that read any output of this task. */
        Set<String> consumerProcesses = Collections.<String>emptySet()
        /** Outputs left in place by the first clean because a consumer may still read them. */
        final Set<String> held = new LinkedHashSet<>()
        /** Outputs that may be published and have not been yet. */
        final Set<String> unpublished = new HashSet<>()
        String dataset
        String sample
        /** The sample subset a per-dataset task was asked for, if it says. */
        List<String> sampleList = Collections.<String>emptyList()
        boolean pinned
        /** The work directory has been cleaned (held outputs aside). */
        boolean deleted
        /** Attempt directories already emptied as failures, which the bulk clean skips. */
        final Set<String> emptiedDirs = new HashSet<>()

        Node(String key, String process, String group) {
            this.key = key
            this.process = process
            this.group = group
        }

        String getWorkDir() { attemptDirs ? attemptDirs.last() : null }

        @Override
        String toString() { "${name ?: key}" }
    }

    private final Map<String,Node> nodes = new LinkedHashMap<>()
    /** work dir -> node, for every attempt */
    private final Map<String,Node> byDir = new HashMap<>()
    /** inputs seen before their producer registered its directory */
    private final Map<String,Set<Node>> orphanInputs = new HashMap<>()
    /** files published before their task completed */
    private final Set<String> earlyPublished = new HashSet<>()
    private final Set<String> terminated = new HashSet<>()
    private final Set<List<String>> pinnedKeys = new HashSet<>()

    private final List<Pattern> barrierPatterns
    private final List<Pattern> singleUsePatterns
    private final List<String> barrierNames
    private final List<String> singleUseNames
    private final Map<String,String> groupCache = new HashMap<>()

    /** Pinning reached a node that had already been deleted: the single-use declaration was wrong. */
    final List<Node> pinnedAfterDelete = []

    TaskGraph(List<String> barrierProcesses, List<String> singleUseProcesses) {
        this.barrierNames = barrierProcesses
        this.singleUseNames = singleUseProcesses
        this.barrierPatterns = barrierProcesses.collect { LineageKeepConfig.wildcard(it) }
        this.singleUsePatterns = singleUseProcesses.collect { LineageKeepConfig.wildcard(it) }
    }

    Collection<Node> getNodes() { nodes.values() }

    Node node(String key) { nodes[key] }

    static String simpleName(String process) {
        final i = process.lastIndexOf(':')
        return i < 0 ? process : process.substring(i + 1)
    }

    /**
     * The group a process belongs to: the barrier or single-use pattern its simple name
     * matches, else the simple name itself.
     */
    String groupOf(String process) {
        groupCache.computeIfAbsent(process) { String p ->
            final name = simpleName(p)
            for( int i = 0; i < barrierPatterns.size(); i++ )
                if( barrierPatterns[i].matcher(name).matches() ) return 'barrier:' + barrierNames[i]
            for( int i = 0; i < singleUsePatterns.size(); i++ )
                if( singleUsePatterns[i].matcher(name).matches() ) return 'single:' + singleUseNames[i]
            return name
        }
    }

    static boolean isBarrier(String group) { group.startsWith('barrier:') }

    static boolean isSingleUse(String group) { group.startsWith('single:') }

    // ---- events ----------------------------------------------------------

    /**
     * A task attempt was created (or found in the cache). `inputDirs` are the work
     * directories its input files live in.
     */
    Node register(String key, String process, String workDir, Collection<String> inputDirs,
                  String name, String hash, String dataset, String sample) {
        Node n = nodes[key]
        if( n == null ) {
            n = new Node(key, process, groupOf(process))
            nodes[key] = n
        }
        n.name = name
        n.hash = hash
        if( dataset != null ) n.dataset = dataset
        if( sample != null ) n.sample = sample
        if( workDir != null && !n.attemptDirs.contains(workDir) ) {
            n.attemptDirs << workDir
            byDir[workDir] = n
            final waiting = orphanInputs.remove(workDir)
            if( waiting )
                waiting.each { Node c -> link(n, c) }
        }
        for( String dir : inputDirs ) {
            final p = byDir[dir]
            if( p != null )
                link(p, n)
            else
                orphanInputs.computeIfAbsent(dir) { new LinkedHashSet<Node>() }.add(n)
        }
        if( rerunKey(n) in pinnedKeys )
            pin(n)
        return n
    }

    private void link(Node producer, Node consumer) {
        if( producer.is(consumer) )
            return
        // pins travel upstream only: consuming a pinned node does not pin the consumer
        producer.consumers << consumer
        consumer.producers << producer
        if( consumer.pinned )
            pin(producer)
    }

    /**
     * @param outputs every output path the task produced -> the processes that read it
     * @param publishable the outputs that may be published
     */
    void succeeded(String key, Map<String,Set<String>> outputs, Collection<String> publishable) {
        final n = nodes[key]
        if( n == null ) return
        n.state = State.SUCCEEDED
        n.outputs = outputs
        final Set<String> procs = new LinkedHashSet<>()
        outputs.values().each { Set<String> v -> procs.addAll(v) }
        n.consumerProcesses = procs
        for( String path : publishable )
            if( !earlyPublished.remove(path) )
                n.unpublished << path
    }

    /**
     * The task's last attempt failed and will not be retried.
     */
    void failed(String key, String exit) {
        final n = nodes[key]
        if( n == null ) return
        n.state = State.FAILED
        n.exit = exit
        if( isBarrier(n.group) )
            return
        final rk = rerunKey(n)
        if( rk != null ) {
            pinnedKeys << rk
            for( Node other : nodes.values() )
                if( rerunKey(other) == rk ) pin(other)
        }
        pin(n)
    }

    /**
     * The work directories to empty now that an attempt of `key` failed. With `retry` it is
     * that attempt alone, when `cleanRetried`; otherwise the task failed for good and it is
     * every attempt of it, when `cleanFailed`. Nothing reads a failed attempt's directory: it
     * emitted no outputs, and a `-resume` re-runs it elsewhere. Pinning is not affected.
     * Each directory is handed out once.
     */
    List<String> failedDirs(String key, String workDir, boolean retry, boolean cleanRetried, boolean cleanFailed) {
        final n = nodes[key]
        if( n == null )
            return Collections.<String>emptyList()
        final List<String> dirs = retry
            ? (cleanRetried && workDir != null ? [workDir] : Collections.<String>emptyList())
            : (cleanFailed ? n.attemptDirs : Collections.<String>emptyList())
        return dirs.findAll { String d -> n.emptiedDirs.add(d) }
    }

    void published(String path) {
        for( Node n : ownersOf(path) ) {
            if( n.unpublished.remove(path) )
                return
        }
        earlyPublished << path
    }

    private List<Node> ownersOf(String path) {
        final result = []
        for( String dir : parentDirs(path) ) {
            final n = byDir[dir]
            if( n != null ) result << n
        }
        return result as List<Node>
    }

    private static List<String> parentDirs(String path) {
        final result = []
        int i = path.length()
        while( (i = path.lastIndexOf('/', i - 1)) > 0 )
            result << path.substring(0, i)
        return result as List<String>
    }

    void terminated(String process) {
        terminated << process
    }

    // ---- pinning -----------------------------------------------------------

    private void pin(Node start) {
        final queue = new ArrayDeque<Node>()
        queue << start
        while( queue ) {
            final n = queue.poll()
            if( n.pinned ) continue
            n.pinned = true
            if( n.deleted ) pinnedAfterDelete << n
            queue.addAll(n.producers)
        }
    }

    static List<String> rerunKey(Node n) {
        n.sample != null && n.dataset != null ? [n.dataset, n.sample] : null
    }

    // ---- deletion ----------------------------------------------------------

    static class Cleanup {
        final Node node
        /** true: clean the work directories except `keep`; false: delete `paths` only */
        final boolean bulk
        final Set<String> paths
        Cleanup(Node node, boolean bulk, Set<String> paths) {
            this.node = node
            this.bulk = bulk
            this.paths = paths
        }
    }

    /**
     * What may be deleted now. A node whose lineage is OK gets one bulk clean, keeping the
     * outputs a consumer may still read; later calls release those as their readers end.
     * With `flowComplete` every process has ended, so publishing is over too.
     */
    List<Cleanup> collectDeletable(boolean flowComplete) {
        final memo = new HashMap<Node,Boolean>()
        final List<Cleanup> result = []
        for( Node n : nodes.values() ) {
            if( n.pinned || n.state != State.SUCCEEDED )
                continue
            if( !n.deleted ) {
                if( !ok(n, memo, flowComplete) )
                    continue
                n.deleted = true
                n.outputs.each { String path, Set<String> procs -> if( !pathFree(n, procs) ) n.held << path }
                result << new Cleanup(n, true, new LinkedHashSet<String>(n.held))
            }
            else if( n.held ) {
                final free = n.held.findAll { String path -> pathFree(n, n.outputs[path]) } as Set<String>
                if( free ) {
                    n.held.removeAll(free)
                    result << new Cleanup(n, false, free)
                }
            }
        }
        return result
    }

    @PackageScope
    boolean ok(Node n, Map<Node,Boolean> memo, boolean flowComplete) {
        final cached = memo[n]
        if( cached != null )
            return cached
        memo[n] = false   // guards against a cycle, which a DAG should never have
        boolean result = n.state == State.SUCCEEDED && !n.pinned
            && (flowComplete || n.unpublished.isEmpty())
        if( result ) {
            for( String g : groupsOf(n.consumerProcesses) ) {
                if( !isBarrier(g) && !satisfied(n, n.consumerProcesses, g) ) { result = false; break }
            }
        }
        if( result ) {
            for( Node c : n.consumers ) {
                if( isBarrier(c.group) ) continue
                if( c.state != State.SUCCEEDED || !ok(c, memo, flowComplete) ) { result = false; break }
            }
        }
        memo[n] = result
        return result
    }

    private Set<String> groupsOf(Collection<String> processes) {
        final result = new LinkedHashSet<String>()
        for( String p : processes ) result << groupOf(p)
        return result
    }

    /**
     * No reader of the group is still to come: either every process of it among `procs` has
     * terminated, or it is single-use and one of its tasks already consumed this node.
     */
    private boolean satisfied(Node n, Collection<String> procs, String group) {
        final inGroup = procs.findAll { String p -> groupOf(p) == group }
        if( inGroup.every { String p -> p in terminated } )
            return true
        return isSingleUse(group) && n.consumers.any { Node c -> c.group == group }
    }

    /** No reader of an output read by `procs` is still to come. */
    private boolean pathFree(Node n, Set<String> procs) {
        if( !procs ) return true
        for( String g : groupsOf(procs) )
            if( !satisfied(n, procs, g) ) return false
        return true
    }

    // ---- reporting ---------------------------------------------------------

    List<Node> finalFailures() {
        nodes.values().findAll { Node n -> n.state == State.FAILED } as List<Node>
    }

    /**
     * Every node upstream of `n`, nearest first.
     */
    static List<Node> ancestors(Node n) {
        final seen = new LinkedHashSet<Node>()
        final queue = new ArrayDeque<Node>(n.producers)
        while( queue ) {
            final p = queue.poll()
            if( seen.add(p) ) queue.addAll(p.producers)
        }
        return seen as List<Node>
    }

    /**
     * The (dataset, sample) pairs to re-run for the final failures. A failure that carries a
     * sample re-runs that sample. One that carries only a dataset (a per-dataset step) re-runs
     * the sample subset it was given, else the samples upstream of it, else the whole dataset
     * (empty sample).
     * Failures of barrier processes, and failures without a dataset, are not re-run items.
     */
    Map<String,Set<String>> rerunItems() {
        final Map<String,Set<String>> result = new TreeMap<>()
        for( Node n : finalFailures() ) {
            if( isBarrier(n.group) || n.dataset == null )
                continue
            final samples = result.computeIfAbsent(n.dataset) { new TreeSet<String>() }
            if( n.sample != null ) {
                samples << n.sample
                continue
            }
            final upstream = ancestors(n).findAll { Node a -> a.dataset == n.dataset && a.sample != null }
            if( n.sampleList )
                samples.addAll(n.sampleList)
            else if( upstream )
                upstream.each { Node a -> samples << a.sample }
            else
                samples << ''
        }
        return result
    }
}
