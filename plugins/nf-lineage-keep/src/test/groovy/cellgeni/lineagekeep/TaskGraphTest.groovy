package cellgeni.lineagekeep

import spock.lang.Specification

/**
 * The pipeline's chain in miniature:
 *
 *   WGET (per run) -> RENAME_RUN (per run) -> RENAME_SAMPLE (per sample) -> STARSOLO (per sample)
 *       -> MAPPINGQC (per dataset) -> COLLATE (barrier, once)
 *   STARSOLO -> RUNREPORT (barrier, once)
 */
class TaskGraphTest extends Specification {

    static final String WGET = 'WF:DL:WGET10X'
    static final String RUN = 'WF:DL:RENAME10XRUN'
    static final String SAMPLE = 'WF:DL:RENAME10XSAMPLE'
    static final String STAR_H = 'WF:STARSOLO10X_HUMAN:REPROCESS10X_STARSOLO10X'
    static final String STAR_M = 'WF:STARSOLO10X_MOUSE:REPROCESS10X_STARSOLO10X'
    static final String QC = 'WF:REPROCESS10X_MAPPINGQC'
    static final String COLLATE = 'REPROCESS10X_COLLATE'
    static final String REPORT = 'REPROCESS10X_RUNREPORT'

    TaskGraph graph

    def setup() {
        graph = new TaskGraph(['REPROCESS10X_COLLATE', 'REPROCESS10X_RUNREPORT'],
                              ['RENAME10XRUN', 'RENAME10XSAMPLE', 'REPROCESS10X_STARSOLO10X', 'REPROCESS10X_MAPPINGQC'])
    }

    /**
     * register + succeed, returning the node's work dir. The task gets one output per
     * consumer process, `<dir>/to_<simple name>`, read by that process alone.
     */
    private String done(String key, String process, List<String> inputs, Set<String> consumers,
                        String dataset = 'GSE1', String sample = 'S1', List<String> publishable = []) {
        final dir = "/w/${key.replaceAll('[^A-Za-z0-9]', '')}"
        graph.register(key, process, dir, inputs, key, 'h' + key, dataset, sample)
        final Map<String,Set<String>> outputs = consumers.collectEntries { String p ->
            [("$dir/to_${TaskGraph.simpleName(p)}".toString()): [p] as Set<String>]
        }
        graph.succeeded(key, outputs, publishable)
        return dir
    }

    /** one sample through STARsolo; QC not yet run */
    private Map<String,String> sampleChain(String s, String run) {
        final w = done("wget/$run", WGET, [], [RUN, COLLATE] as Set, 'GSE1', s)
        final r = done("run/$run", RUN, [w], [SAMPLE, COLLATE] as Set, 'GSE1', s)
        final m = done("sample/$s", SAMPLE, [r], [STAR_H, STAR_M] as Set, 'GSE1', s)
        final a = done("star/$s", STAR_H, [m], [QC, REPORT] as Set, 'GSE1', s)
        return [wget: w, run: r, sample: m, star: a]
    }

    /** keys of the nodes whose work dir gets its first (bulk) clean now */
    private Set<String> deletableKeys(boolean complete = false) {
        graph.collectDeletable(complete).findAll { it.bulk }.collect { it.node.key } as Set
    }

    def 'groups processes by matching pattern, so species aliases share one group'() {
        expect:
        graph.groupOf(STAR_H) == graph.groupOf(STAR_M)
        TaskGraph.isSingleUse(graph.groupOf(STAR_H))
        TaskGraph.isBarrier(graph.groupOf(REPORT))
        graph.groupOf(WGET) == 'WGET10X'
    }

    def 'nothing is deleted while the lineage is still running'() {
        given:
        sampleChain('S1', 'R1')

        expect: 'STARsolo succeeded but its QC has not run'
        deletableKeys().isEmpty()
    }

    def 'a fully successful lineage is deleted up to the barrier'() {
        given:
        final c = sampleChain('S1', 'R1')
        done('qc/GSE1', QC, [c.star], [COLLATE] as Set, 'GSE1', null)

        when:
        final first = graph.collectDeletable(false)

        then: 'every work dir is cleaned, but what the barriers will read is held'
        first.every { it.bulk }
        first.collect { it.node.key } as Set == ['wget/R1', 'run/R1', 'sample/S1', 'star/S1', 'qc/GSE1'] as Set
        first.find { it.node.key == 'star/S1' }.paths == ['/w/starS1/to_REPROCESS10X_RUNREPORT'] as Set
        first.find { it.node.key == 'wget/R1' }.paths == ['/w/wgetR1/to_REPROCESS10X_COLLATE'] as Set
        first.find { it.node.key == 'sample/S1' }.paths.isEmpty()

        when:
        graph.terminated(COLLATE)
        graph.terminated(REPORT)
        final second = graph.collectDeletable(false)

        then: 'the held outputs go once the barriers have ended'
        second.every { !it.bulk }
        second.collectMany { it.paths } as Set == [
            '/w/wgetR1/to_REPROCESS10X_COLLATE', '/w/runR1/to_REPROCESS10X_COLLATE',
            '/w/starS1/to_REPROCESS10X_RUNREPORT', '/w/qcGSE1/to_REPROCESS10X_COLLATE'] as Set
        graph.collectDeletable(true).isEmpty()
    }

    def 'a final failure keeps the failed task and everything upstream of it'() {
        given:
        final w = done('wget/R1', WGET, [], [RUN, COLLATE] as Set)
        final r = done('run/R1', RUN, [w], [SAMPLE, COLLATE] as Set)
        final m = done('sample/S1', SAMPLE, [r], [STAR_H, STAR_M] as Set)
        graph.register('star/S1', STAR_H, '/w/star1', [m], 'star/S1', 'h', 'GSE1', 'S1')

        when:
        graph.failed('star/S1', '1')
        [STAR_H, STAR_M, QC, COLLATE, REPORT, RUN, SAMPLE, WGET].each { graph.terminated(it) }

        then:
        deletableKeys(true).isEmpty()
        graph.nodes.findAll { it.pinned }*.key as Set == ['wget/R1', 'run/R1', 'sample/S1', 'star/S1'] as Set
        graph.rerunItems() == [GSE1: ['S1'] as Set]
    }

    def 'a sibling sample of the failed one is still cleaned'() {
        given:
        final ok = sampleChain('S2', 'R2')
        final w = done('wget/R1', WGET, [], [RUN, COLLATE] as Set, 'GSE1', 'S1')
        graph.register('run/R1', RUN, '/w/run1', [w], 'run/R1', 'h', 'GSE1', 'S1')
        graph.failed('run/R1', '1')
        done('qc/GSE1', QC, [ok.star], [COLLATE] as Set, 'GSE1', null)
        [RUN, SAMPLE].each { graph.terminated(it) }

        expect:
        deletableKeys() == ['wget/R2', 'run/R2', 'sample/S2', 'star/S2', 'qc/GSE1'] as Set
        graph.node('wget/R1').pinned
    }

    def 'a failed run pins the sibling runs of its sample, even when the sample went ahead without it'() {
        given: 'R1 fails; R2 of the same sample carries on and aligns on partial data'
        final w1 = done('wget/R1', WGET, [], [RUN, COLLATE] as Set, 'GSE1', 'S1')
        graph.register('run/R1', RUN, '/w/run1', [w1], 'run/R1', 'h', 'GSE1', 'S1')
        graph.failed('run/R1', '1')
        final c = sampleChain('S1', 'R2')
        done('qc/GSE1', QC, [c.star], [COLLATE] as Set, 'GSE1', null)
        [COLLATE, REPORT].each { graph.terminated(it) }

        expect:
        deletableKeys(true) == ['qc/GSE1'] as Set
        graph.node('wget/R2').pinned
        graph.node('star/S1').pinned
    }

    def 'a retried attempt releases nothing and the retry is tracked in the same node'() {
        given:
        final w = done('wget/R1', WGET, [], [RUN] as Set)
        graph.register('run/R1', RUN, '/w/run1a', [w], 'run/R1', 'h1', 'GSE1', 'S1')

        when: 'attempt 1 fails with RETRY: the observer records nothing; attempt 2 registers'
        graph.register('run/R1', RUN, '/w/run1b', [w], 'run/R1', 'h2', 'GSE1', 'S1')

        then:
        graph.node('run/R1').attemptDirs == ['/w/run1a', '/w/run1b']
        deletableKeys().isEmpty()
    }

    def 'failed attempts are handed out for emptying per switch, once each, and pinning is unchanged'() {
        given:
        final w = done('wget/R1', WGET, [], [RUN] as Set)
        graph.register('run/R1', RUN, '/w/run1a', [w], 'run/R1', 'h1', 'GSE1', 'S1')

        expect: 'a retried attempt goes only with cleanRetried'
        graph.failedDirs('run/R1', '/w/run1a', true, false, true).isEmpty()
        graph.failedDirs('run/R1', '/w/run1a', true, true, false) == ['/w/run1a']
        graph.failedDirs('run/R1', '/w/run1a', true, true, false).isEmpty()

        when: 'attempt 2 fails for good'
        graph.register('run/R1', RUN, '/w/run1b', [w], 'run/R1', 'h2', 'GSE1', 'S1')
        graph.failed('run/R1', '1')

        then: 'cleanFailed empties every attempt not emptied yet; the task and its producer stay pinned'
        graph.failedDirs('run/R1', '/w/run1b', false, true, false).isEmpty()
        graph.failedDirs('run/R1', '/w/run1b', false, false, true) == ['/w/run1b']
        graph.node('run/R1').pinned
        graph.node('wget/R1').pinned
        deletableKeys(true).isEmpty()
    }

    def 'with cleanFailed alone, a final failure empties its earlier retried attempts too'() {
        given:
        graph.register('wget/R1', WGET, '/w/w1a', [], 'wget/R1', 'h1', 'GSE1', 'S1')
        graph.register('wget/R1', WGET, '/w/w1b', [], 'wget/R1', 'h2', 'GSE1', 'S1')
        graph.failed('wget/R1', '1')

        expect:
        graph.failedDirs('wget/R1', '/w/w1b', false, false, true) == ['/w/w1a', '/w/w1b']
    }

    def 'the bulk clean of a task whose retry succeeded skips the attempt already emptied'() {
        given:
        graph.register('wget/R1', WGET, '/w/w1a', [], 'wget/R1', 'h1', 'GSE1', 'S1')
        graph.failedDirs('wget/R1', '/w/w1a', true, true, false)
        graph.register('wget/R1', WGET, '/w/w1b', [], 'wget/R1', 'h2', 'GSE1', 'S1')
        graph.succeeded('wget/R1', [:], [])

        when:
        final c = graph.collectDeletable(true)

        then: 'the observer cleans attemptDirs minus emptiedDirs'
        c*.node*.key == ['wget/R1']
        c[0].node.attemptDirs - c[0].node.emptiedDirs == ['/w/w1b']
    }

    def 'a non-single-use consumer process must terminate first'() {
        given:
        final w = done('wget/R1', WGET, [], ['WF:OTHER'] as Set)
        done('other/1', 'WF:OTHER', [w], [] as Set)

        expect: 'the leaf goes at once; its producer waits, since OTHER may still start another task on it'
        deletableKeys() == ['other/1'] as Set

        when:
        graph.terminated('WF:OTHER')

        then:
        deletableKeys() == ['wget/R1'] as Set
    }

    def 'unpublished outputs hold a task until published, or until the run completes'() {
        given:
        done('wget/R1', WGET, [], [] as Set, 'GSE1', 'S1', ['/w/wgetR1/R1.fastq.gz'])

        expect:
        deletableKeys().isEmpty()

        when:
        graph.published('/w/wgetR1/R1.fastq.gz')

        then:
        deletableKeys() == ['wget/R1'] as Set
    }

    def 'a file published before its task completed still counts'() {
        given:
        graph.register('wget/R1', WGET, '/w/wgetR1', [], 'wget/R1', 'h', 'GSE1', 'S1')
        graph.published('/w/wgetR1/R1.fastq.gz')
        graph.succeeded('wget/R1', [:], ['/w/wgetR1/R1.fastq.gz'])

        expect:
        deletableKeys() == ['wget/R1'] as Set
    }

    def 'a barrier failure pins nothing and is not a re-run item'() {
        given:
        final c = sampleChain('S1', 'R1')
        done('qc/GSE1', QC, [c.star], [COLLATE] as Set, 'GSE1', null)
        graph.register('report', REPORT, '/w/report', [c.star], 'report', 'h', null, null)
        graph.failed('report', '1')
        [COLLATE, REPORT].each { graph.terminated(it) }

        expect:
        deletableKeys(true).size() == 5
        graph.collectDeletable(true).isEmpty()
        graph.rerunItems().isEmpty()
    }

    def 'a per-dataset failure re-runs its subset, else its upstream samples, else the whole dataset'() {
        given:
        graph.register('fetch/GSE2', 'WF:FETCH10XMETA', '/w/f2', [], 'fetch/GSE2', 'h', 'GSE2', null)
        graph.node('fetch/GSE2').sampleList = ['A', 'B']
        graph.failed('fetch/GSE2', '1')

        graph.register('fetch/GSE3', 'WF:FETCH10XMETA', '/w/f3', [], 'fetch/GSE3', 'h', 'GSE3', null)
        graph.failed('fetch/GSE3', '1')

        final c = sampleChain('S1', 'R1')
        graph.register('qc/GSE1', QC, '/w/qc1', [c.star], 'qc/GSE1', 'h', 'GSE1', null)
        graph.failed('qc/GSE1', '1')

        expect:
        graph.rerunItems() == [GSE1: ['S1'] as Set, GSE2: ['A', 'B'] as Set, GSE3: [''] as Set]
        graph.node('wget/R1').pinned
    }

    def 'an input seen before its producer registered is linked once it does'() {
        given:
        graph.register('run/R1', RUN, '/w/run1', ['/w/wget1'], 'run/R1', 'h', 'GSE1', 'S1')
        graph.register('wget/R1', WGET, '/w/wget1', [], 'wget/R1', 'h', 'GSE1', 'S1')

        expect:
        graph.node('run/R1').producers*.key == ['wget/R1']
    }

    def 'parses the task directory of a path under the work root'() {
        expect:
        LineageKeepObserver.taskDirOf('/x/work/ab/cdef/sub/f.gz', '/x/work') == '/x/work/ab/cdef'
        LineageKeepObserver.taskDirOf('/x/work/ab/cdef', '/x/work') == '/x/work/ab/cdef'
        LineageKeepObserver.taskDirOf('/x/work/ab', '/x/work') == null
        LineageKeepObserver.taskDirOf('/nfs/ref/index', '/x/work') == null
        LineageKeepObserver.taskDirOf('/x/workdir/ab/cd', '/x/work') == null
    }
}
