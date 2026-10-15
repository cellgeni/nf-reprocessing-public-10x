import groovy.util.logging.Slf4j
import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.ConcurrentHashMap
import java.util.zip.GZIPInputStream

/**
 * STARsolo memory, run limit and queue predicted from the sample's input. plan() runs in
 * subworkflows/local/starsolo10x and its result travels in `meta.resources`, which the
 * `memory`, `time` and `queue` closures in configs/starsolo10x.config read.
 *
 * Peak memory and runtime both grow with STAR's input read count. Nothing upstream knows that
 * count, so it is estimated from the R2 files: the head of each is streamed through gunzip, which gives
 * compressed bytes per read and mean read length, and the file size divided by bytes-per-read
 * gives the read count. Only R2 is read because STAR's "Number of input reads" and "Average input
 * read length" are the cDNA read's.
 *
 * The coefficients and the sampling depth come from notebooks/resource_prediction.ipynb.
 * Refit there and paste the printed constants below. Any failure (no R2, unreadable file) falls
 * back to the old fixed memory ladder and 48 h on `long` rather than failing the task.
 */
@Slf4j
class StarsoloResources {

    // Fitted in notebooks/resource_prediction.ipynb (§10) on batches 9-21, n=11,998:
    // q0.99 quantile regressions of peak RSS and realtime on input reads. On hold-out batches
    // 22-24 this policy gave no memlimit or runlimit kills and sent ~95% of samples to `normal`.
    static final double BASE_GB        = 37.7
    static final double GB_PER_GREAD   = 111.2
    static final double MARGIN         = 1.1
    static final int    FLOOR_GB       = 48
    static final int    CAP_GB         = 1024
    static final double RETRY_STEP     = 1.5
    static final int    SAMPLE_RECORDS = 50_000
    static final double BASE_H         = 0.13
    static final double H_PER_GREAD    = 7.13
    static final double TIME_MARGIN    = 2.0

    // LSF queues: same priority and a shared fairshare, but only `normal` backfills, and only a job
    // with an explicit run limit short enough to finish before a reservation can be backfilled.
    static final int    NORMAL_LIMIT_H = 12
    static final int    LONG_LIMIT_H   = 48
    // configs/starsolo10x.config's maxRetries; errorStrategy stops at 3, the rest is headroom
    static final int    MAX_ATTEMPTS   = 5

    // Mouse has no history yet (every run so far was human) and its 2020A index is the same
    // size class, so it uses the human fit until there is data to refit it.
    static final Map<String, Double> REF_BASE_GB = [human: BASE_GB, mouse: BASE_GB]

    // Sample each file once per run, whatever calls in.
    private static final Map<String, Map> CACHE = new ConcurrentHashMap<>()

    /** Counts the compressed bytes GZIPInputStream pulls from the file. */
    private static class CountingInputStream extends FilterInputStream {
        long count = 0
        CountingInputStream(InputStream in) { super(in) }
        @Override int read() {
            int b = super.read()
            if (b >= 0) count++
            return b
        }
        @Override int read(byte[] buf, int off, int len) {
            int n = super.read(buf, off, len)
            if (n > 0) count += n
            return n
        }
        @Override long skip(long n) {
            long s = super.skip(n)
            count += s
            return s
        }
    }

    /** Unwrap a TaskPath (a staged input) to the real file; plain paths pass through. */
    private static Path realPath(Object p) {
        def path = p.respondsTo('getTarget') ? p.getTarget() : p
        return path as Path
    }

    /**
     * Estimate one gzipped FASTQ's read count and mean read length from its first `records` reads.
     * A file shorter than that is counted exactly.
     */
    static Map sampleFastq(Path gz, int records = SAMPLE_RECORDS) {
        def key = "${gz.toUri()}:${records}".toString()
        def hit = CACHE.get(key)
        if (hit != null)
            return hit

        long size = Files.size(gz)
        long recs = 0, seqLen = 0, lines = 0
        def counter = new CountingInputStream(Files.newInputStream(gz))
        boolean eof = false
        try {
            def reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(counter, 8192), 'US-ASCII'))
            String line
            while (recs < records) {
                line = reader.readLine()
                if (line == null) {
                    eof = true
                    break
                }
                if (lines % 4 == 1)
                    seqLen += line.length()
                lines++
                if (lines % 4 == 0)
                    recs++
            }
        }
        finally {
            counter.close()
        }
        if (recs == 0)
            throw new IllegalStateException("no FASTQ records in ${gz}")
        double reads = eof ? recs : size / (counter.count / (double) recs)
        def result = [reads: reads, readLength: seqLen / (double) recs, exact: eof]
        CACHE.put(key, result)
        return result
    }

    /** Estimated input reads and mean read length of a sample, summed over its R2 files. */
    static Map estimateReads(Object fastqs, int records = SAMPLE_RECORDS) {
        def files = (fastqs instanceof Collection ? fastqs : [fastqs]) as List
        def r2 = files.findAll { it.getFileName().toString() =~ /_R2_/ }
        if (!r2)
            throw new IllegalArgumentException("no R2 FASTQ among ${files*.getFileName()}")
        double reads = 0, lengthSum = 0
        for (f in r2) {
            def s = sampleFastq(realPath(f), records)
            reads += s.reads
            lengthSum += s.reads * s.readLength
        }
        return [reads: reads, readLength: reads ? lengthSum / reads : 0d, files: r2.size()]
    }

    /** Memory request in GB for one attempt, given an input read count. */
    static long requestGb(double reads, String ref, int attempt) {
        double base = REF_BASE_GB.getOrDefault(ref, BASE_GB)
        double first = Math.max(FLOOR_GB, (base + GB_PER_GREAD * reads / 1e9) * MARGIN)
        double req = roundUp16(roundUp16(first) * Math.pow(RETRY_STEP, attempt - 1))
        return (long) Math.min(CAP_GB, req)
    }

    /** Run limit in whole hours for one attempt, given an input read count; doubles on each retry. */
    static long limitHours(double reads, int attempt) {
        double first = Math.ceil(Math.max(1d, (BASE_H + H_PER_GREAD * reads / 1e9) * TIME_MARGIN))
        return (long) Math.min(LONG_LIMIT_H, first * Math.pow(2, attempt - 1))
    }

    /**
     * Per-attempt resources for one STARsolo sample, for `meta.resources`:
     * `[reads:, mem_gb: [...], hours: [...], queue: [...]]`, one list entry per attempt.
     * configs/starsolo10x.config indexes the lists with `task.attempt - 1`; a config closure
     * cannot see `lib/` classes, so the work is done here, in the workflow, once per sample.
     * Returns null when the estimate fails, which the config reads as "use the fixed ladder".
     */
    static Map plan(Object fastqs, String ref) {
        try {
            def est = estimateReads(fastqs)
            double reads = est.reads as double
            def attempts = 1..MAX_ATTEMPTS
            def hours = attempts.collect { limitHours(reads, it) }
            log.debug "STARsolo resources: ${String.format('%.3g', reads / 1e9)}e9 reads " +
                      "(${est.files} R2 files, mean length ${Math.round(est.readLength as double)})"
            return [
                reads : Math.round(reads),
                mem_gb: attempts.collect { requestGb(reads, ref, it) },
                hours : hours,
                queue : hours.collect { it <= NORMAL_LIMIT_H ? 'normal' : 'long' },
            ]
        }
        catch (Exception e) {
            log.warn "STARsolo resource estimate failed (${e.message}); using the fixed ladder on long"
            return null
        }
    }

    private static double roundUp16(double gb) {
        return Math.ceil(gb / 16d) * 16d
    }
}
