package cellgeni.lineagekeep

import java.nio.file.FileSystems
import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes
import java.time.Instant

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Empties work directories and records each one in `deleted.tsv`.
 *
 * Sizes are counted while walking, so a dry run reports the same numbers a real run would.
 * `bytes` is the apparent size; `freed_bytes` counts only files with a single link, because
 * a work file that was published with `mode 'link'` (or that a rename step hard-linked) keeps
 * its blocks for as long as the other link exists.
 */
@Slf4j
@CompileStatic
class WorkDirCleaner implements Closeable {

    static final String HEADER = 'time\taction\tprocess\tname\thash\tworkdir\tfiles\tbytes\tfreed_bytes'

    private final boolean dryRun
    private final boolean keepLogs
    private final Writer log0

    long totalBytes
    long totalFreed
    int dirs

    WorkDirCleaner(Path tsv, boolean dryRun, boolean keepLogs) {
        this.dryRun = dryRun
        this.keepLogs = keepLogs
        if( tsv != null ) {
            Files.createDirectories(tsv.parent)
            final exists = Files.exists(tsv)
            this.log0 = Files.newBufferedWriter(tsv, java.nio.file.StandardOpenOption.CREATE, java.nio.file.StandardOpenOption.APPEND)
            if( !exists ) { log0.write(HEADER); log0.write('\n'); log0.flush() }
        }
    }

    static class Stats {
        int files
        long bytes
        long freed
    }

    /**
     * Delete everything in `dir` except `keep` (absolute paths) and, with `keepLogs`, the
     * top-level `.command.*` and `.exitcode` files. The directory itself stays, so the
     * trace's workdir column still points somewhere.
     */
    Stats clean(Path dir, Set<String> keep, String process, String name, String hash) {
        return empty(dir, keep ?: Collections.<String>emptySet(), -1, dryRun ? 'dry-run' : 'delete', process, name, hash)
    }

    /**
     * Empty the work directory of a failed attempt. Besides the logs, regular files smaller
     * than `keepBelow` bytes and symlinks survive: they cost nothing and are what the failure
     * is debugged from. `kind` (`retried` or `failed`) is the action recorded.
     */
    Stats cleanFailed(Path dir, long keepBelow, String kind, String process, String name, String hash) {
        return empty(dir, Collections.<String>emptySet(), keepBelow, dryRun ? "dry-run-${kind}".toString() : kind, process, name, hash)
    }

    /** `keepBelow` < 0: no size rule, symlinks go too */
    private Stats empty(Path dir, Set<String> keep, long keepBelow, String action, String process, String name, String hash) {
        final stats = new Stats()
        if( dir == null || !Files.isDirectory(dir, LinkOption.NOFOLLOW_LINKS) )
            return stats
        try ( final entries = Files.newDirectoryStream(dir) ) {
            for( Path entry : entries ) {
                if( keepLogs && isLog(entry.fileName.toString()) )
                    continue
                removeTree(entry, keep, stats, keepBelow)
            }
        }
        catch( Exception e ) {
            log.warn "[nf-lineage-keep] unable to clean ${dir} -- ${e.message ?: e}"
        }
        record(action, process, name, hash, dir, stats)
        return stats
    }

    /**
     * Delete the given paths (outputs released after the first clean).
     */
    Stats remove(Collection<String> paths, String process, String name, String hash) {
        final stats = new Stats()
        for( String p : paths ) {
            try {
                final path = Path.of(p)
                if( Files.exists(path, LinkOption.NOFOLLOW_LINKS) )
                    removeTree(path, Collections.<String>emptySet(), stats)
            }
            catch( Exception e ) {
                log.warn "[nf-lineage-keep] unable to delete ${p} -- ${e.message ?: e}"
            }
        }
        if( paths )
            record(dryRun ? 'dry-run-held' : 'delete-held', process, name, hash, Path.of(paths.first()).parent, stats)
        return stats
    }

    /**
     * Remove whatever under `dir` matches one of the globs (relative to `dir`). Used for
     * scratch such as STAR's `_STARtmp`, at any outcome of the task.
     */
    Stats purge(Path dir, List<String> globs, String process, String name, String hash) {
        final stats = new Stats()
        if( dir == null || !globs || !Files.isDirectory(dir, LinkOption.NOFOLLOW_LINKS) )
            return stats
        final List<PathMatcher> matchers = globs.collect { String g -> FileSystems.default.getPathMatcher('glob:' + g) }
        final List<Path> hits = []
        try {
            Files.walkFileTree(dir, EnumSet.noneOf(java.nio.file.FileVisitOption), 6, new SimpleFileVisitor<Path>() {
                @Override
                FileVisitResult preVisitDirectory(Path d, BasicFileAttributes attrs) {
                    if( d != dir && matchers.any { PathMatcher m -> m.matches(dir.relativize(d)) } ) {
                        hits << d
                        return FileVisitResult.SKIP_SUBTREE
                    }
                    return FileVisitResult.CONTINUE
                }
                @Override
                FileVisitResult visitFile(Path f, BasicFileAttributes attrs) {
                    if( matchers.any { PathMatcher m -> m.matches(dir.relativize(f)) } )
                        hits << f
                    return FileVisitResult.CONTINUE
                }
                @Override
                FileVisitResult visitFileFailed(Path f, IOException e) {
                    return FileVisitResult.CONTINUE
                }
            })
            for( Path p : hits )
                removeTree(p, Collections.<String>emptySet(), stats)
        }
        catch( Exception e ) {
            log.warn "[nf-lineage-keep] unable to purge scratch in ${dir} -- ${e.message ?: e}"
        }
        if( hits )
            record(dryRun ? 'scratch-dry-run' : 'scratch', process, name, hash, dir, stats)
        return stats
    }

    static boolean isLog(String fileName) {
        fileName.startsWith('.command.') || fileName == '.exitcode'
    }

    private void removeTree(Path root, Set<String> keep, Stats stats, long keepBelow = -1) {
        Files.walkFileTree(root, new SimpleFileVisitor<Path>() {
            @Override
            FileVisitResult preVisitDirectory(Path d, BasicFileAttributes attrs) {
                return d.toString() in keep ? FileVisitResult.SKIP_SUBTREE : FileVisitResult.CONTINUE
            }
            @Override
            FileVisitResult visitFile(Path f, BasicFileAttributes attrs) {
                if( f.toString() in keep )
                    return FileVisitResult.CONTINUE
                if( keepBelow >= 0 && (attrs.isSymbolicLink() || (attrs.isRegularFile() && attrs.size() < keepBelow)) )
                    return FileVisitResult.CONTINUE
                stats.files++
                if( attrs.isRegularFile() ) {
                    stats.bytes += attrs.size()
                    if( nlink(f) <= 1 ) stats.freed += attrs.size()
                }
                if( !dryRun ) deleteQuietly(f)
                return FileVisitResult.CONTINUE
            }
            @Override
            FileVisitResult visitFileFailed(Path f, IOException e) {
                if( !dryRun ) deleteQuietly(f)
                return FileVisitResult.CONTINUE
            }
            @Override
            FileVisitResult postVisitDirectory(Path d, IOException e) {
                // a directory still holding a kept path is not empty, and stays
                if( !dryRun && !keep.any { String k -> k.startsWith(d.toString() + '/') } )
                    deleteQuietly(d)
                return FileVisitResult.CONTINUE
            }
        })
    }

    private static int nlink(Path f) {
        try {
            return Files.getAttribute(f, 'unix:nlink', LinkOption.NOFOLLOW_LINKS) as int
        }
        catch( Exception ignored ) {
            return 1
        }
    }

    private static void deleteQuietly(Path p) {
        try {
            Files.deleteIfExists(p)
        }
        catch( NoSuchFileException ignored ) { }
        catch( Exception e ) {
            log.debug "[nf-lineage-keep] unable to delete ${p} -- ${e.message ?: e}"
        }
    }

    private synchronized void record(String action, String process, String name, String hash, Path dir, Stats stats) {
        totalBytes += stats.bytes
        totalFreed += stats.freed
        dirs++
        if( log0 == null ) return
        final row = [Instant.now().toString(), action, process, name, hash, dir.toString(),
                     stats.files.toString(), stats.bytes.toString(), stats.freed.toString()]
        log0.write(row.collect { String s -> (s ?: '').replace('\t', ' ') }.join('\t'))
        log0.write('\n')
        log0.flush()
    }

    @Override
    void close() {
        log0?.close()
    }
}
