package cellgeni.lineagekeep

import java.util.regex.Pattern

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration
import nextflow.util.MemoryUnit

@ScopeName('lineageKeep')
@Description('''
    The `lineageKeep` scope configures the `nf-lineage-keep` plugin, which deletes the work
    directories of tasks whose whole downstream lineage succeeded, keeps every task upstream of
    a failure, and writes a table of the failed items for a re-run.
''')
@CompileStatic
class LineageKeepConfig implements ConfigScope {

    @ConfigOption
    @Description('Enable the plugin (default: `false`).')
    final boolean enabled

    @ConfigOption
    @Description('''
        Log what would be deleted to `deleted.tsv` without deleting anything (default: `true`).
    ''')
    final boolean dryRun

    @ConfigOption
    @Description('''
        Keep `.command.*` and `.exitcode` in a cleaned work directory, so the trace and the
        task logs still resolve (default: `true`).
    ''')
    final boolean keepLogs

    @ConfigOption
    @Description('''
        Keep the work directories of failed attempts whose retry succeeded when the task is
        cleaned (default: `false`). Has no effect on attempts `cleanRetriedAttempts` emptied.
    ''')
    final boolean keepFailedAttempts

    @ConfigOption
    @Description('''
        Empty an attempt that failed and will be retried as soon as it ends, instead of with
        its task once the lineage is done (default: `false`). The directory itself stays: a
        `-resume` only reaches attempt n+1's cache through attempt n's directory.
    ''')
    final boolean cleanRetriedAttempts

    @ConfigOption
    @Description('''
        Empty every attempt of a task that failed for good -- ignored, aborted, or a barrier
        (default: `false`). A `-resume` re-runs such a task in a new directory, so nothing in
        them is reused. The task, its ancestors and its sample stay pinned; only the failed
        task's own directories are emptied, and the directories themselves stay.
    ''')
    final boolean cleanFailedTasks

    @ConfigOption
    @Description('''
        What survives `cleanRetriedAttempts` and `cleanFailedTasks`, besides the logs:
        regular files smaller than this, and symlinks (default: `1 MB`). Keeps the small
        reports a failure is debugged from, such as STAR's `Log.out` or a `chemistry.json`.
    ''')
    final MemoryUnit failedKeepBelow

    @ConfigOption
    @Description('''
        Processes that collect everything at the end of the run (reports, run-level tables).
        Their completion is waited for before deleting what they read, but their success is not
        required for an upstream lineage to count as successful, and their failure pins nothing.
        Matched against the process simple name; `*` is a wildcard.
    ''')
    final List<String> barrierProcesses

    @ConfigOption
    @Description('''
        Processes in which each upstream file is consumed by exactly one task. Once a task of
        such a process has consumed a file, no other task of it is waited for. Without this the
        process has to terminate first. Processes matching one pattern form one group, so
        `STARSOLO*` covers the per-species aliases together. `*` is a wildcard.
    ''')
    final List<String> singleUseProcesses

    @ConfigOption
    @Description('''
        Scratch left in a work directory, removed as soon as any attempt of the task ends,
        whatever its outcome, as `PROCESS:glob` (e.g. `'*STARSOLO10X:**/_STARtmp'`).
    ''')
    final List<String> scratch

    @ConfigOption
    @Description('''
        Name of the task input holding the `meta` map that the re-run table is built from
        (default: `meta`).
    ''')
    final String metaVariable

    @ConfigOption
    @Description('Key of the dataset id in `meta` (default: `dataset_id`).')
    final String datasetKey

    @ConfigOption
    @Description('Key of the sample id in `meta` (default: `sample_id`).')
    final String sampleKey

    @ConfigOption
    @Description('Key of the item id in `meta` (default: `id`).')
    final String idKey

    @ConfigOption
    @Description('''
        Task input holding the sample subset a per-dataset task was asked for, as a list or a
        comma-separated string (default: `sample_ids`). Used for a failure with no sample and
        nothing upstream, so its dataset is re-run with the same subset instead of whole.
    ''')
    final String sampleListVariable

    @ConfigOption
    @Description('Directory for `deleted.tsv`, `failed_tasks.tsv`, `rerun.tsv` (default: `<outputDir>/cleanup`).')
    final String outputDir

    @ConfigOption
    @Description('How often pending events are processed and deletions considered (default: `30s`).')
    final Duration interval

    /* required by extension point -- do not remove */
    LineageKeepConfig() {}

    LineageKeepConfig(Map opts) {
        enabled            = opts.enabled as boolean
        dryRun             = opts.dryRun == null ? true : opts.dryRun as boolean
        keepLogs           = opts.keepLogs == null ? true : opts.keepLogs as boolean
        keepFailedAttempts = opts.keepFailedAttempts as boolean
        cleanRetriedAttempts = opts.cleanRetriedAttempts as boolean
        cleanFailedTasks   = opts.cleanFailedTasks as boolean
        failedKeepBelow    = opts.failedKeepBelow != null ? MemoryUnit.of(opts.failedKeepBelow.toString()) : MemoryUnit.of('1 MB')
        barrierProcesses   = stringList(opts.barrierProcesses)
        singleUseProcesses = stringList(opts.singleUseProcesses)
        scratch            = stringList(opts.scratch)
        metaVariable       = opts.metaVariable ?: 'meta'
        datasetKey         = opts.datasetKey ?: 'dataset_id'
        sampleKey          = opts.sampleKey ?: 'sample_id'
        idKey              = opts.idKey ?: 'id'
        sampleListVariable = opts.sampleListVariable ?: 'sample_ids'
        outputDir          = opts.outputDir as String
        interval           = opts.interval ? Duration.of(opts.interval.toString()) : Duration.of('30s')
    }

    private static List<String> stringList(Object value) {
        if( value == null )
            return Collections.<String>emptyList()
        if( value instanceof Collection )
            return ((Collection)value).collect { it.toString() }
        return value.toString().tokenize(',')*.trim()
    }

    /**
     * Convert a pattern with `*` wildcards into an anchored regex.
     */
    static Pattern wildcard(String glob) {
        final parts = glob.split(/\*/, -1).collect { String s -> Pattern.quote(s) }
        return Pattern.compile(parts.join('.*'))
    }
}
