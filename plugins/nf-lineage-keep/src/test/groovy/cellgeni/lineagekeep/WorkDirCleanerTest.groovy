package cellgeni.lineagekeep

import java.nio.file.Files
import java.nio.file.Path

import spock.lang.Specification
import spock.lang.TempDir

class WorkDirCleanerTest extends Specification {

    @TempDir
    Path tmp

    private Path workDir() {
        final d = Files.createDirectories(tmp.resolve('ab/cdef'))
        Files.writeString(d.resolve('.command.sh'), 'echo hi')
        Files.writeString(d.resolve('.command.err'), 'boom')
        Files.writeString(d.resolve('.exitcode'), '0')
        Files.writeString(d.resolve('out.fastq.gz'), 'x' * 100)
        Files.createDirectories(d.resolve('S1/_STARtmp/sub'))
        Files.writeString(d.resolve('S1/_STARtmp/sub/big'), 'y' * 1000)
        Files.writeString(d.resolve('S1/Log.final.out'), 'z' * 10)
        // a staged input: a symlink to a file elsewhere, which must survive
        final ext = Files.writeString(tmp.resolve('upstream.fastq.gz'), 'u' * 50)
        Files.createSymbolicLink(d.resolve('in.fastq.gz'), ext)
        // a hard link, as `publish mode 'link'` leaves behind
        Files.createLink(tmp.resolve('published.fastq.gz'), d.resolve('out.fastq.gz'))
        return d
    }

    def 'cleans everything but the logs, counting hard-linked bytes as not freed'() {
        given:
        final d = workDir()
        final tsv = tmp.resolve('log/deleted.tsv')
        final cleaner = new WorkDirCleaner(tsv, false, true)

        when:
        final s = cleaner.clean(d, [] as Set, 'P', 'P (1)', 'ab/cdef')
        cleaner.close()

        then:
        Files.list(d).collect { it.fileName.toString() }.sort() == ['.command.err', '.command.sh', '.exitcode']
        Files.exists(tmp.resolve('upstream.fastq.gz'))          // symlink target untouched
        Files.exists(tmp.resolve('published.fastq.gz'))         // other hard link untouched
        s.bytes == 1110
        s.freed == 1010
        Files.readAllLines(tsv).size() == 2
    }

    def 'a dry run deletes nothing but reports the same numbers'() {
        given:
        final d = workDir()
        final cleaner = new WorkDirCleaner(null, true, true)

        when:
        final s = cleaner.clean(d, [] as Set, 'P', 'P (1)', 'ab/cdef')

        then:
        Files.exists(d.resolve('out.fastq.gz'))
        Files.exists(d.resolve('S1/_STARtmp/sub/big'))
        s.bytes == 1110
    }

    def 'purges scratch at any depth and leaves the rest'() {
        given:
        final d = workDir()
        final cleaner = new WorkDirCleaner(null, false, true)

        when:
        final s = cleaner.purge(d, ['**/_STARtmp'], 'P', 'P (1)', 'ab/cdef')

        then:
        !Files.exists(d.resolve('S1/_STARtmp'))
        Files.exists(d.resolve('S1/Log.final.out'))
        Files.exists(d.resolve('out.fastq.gz'))
        s.bytes == 1000
    }

    def 'keepLogs false removes the logs too'() {
        given:
        final d = workDir()
        final cleaner = new WorkDirCleaner(null, false, false)

        when:
        cleaner.clean(d, [] as Set, 'P', 'P (1)', 'ab/cdef')

        then:
        Files.list(d).count() == 0
    }

    def 'held outputs survive the clean, and their parent directories with them, until removed'() {
        given:
        final d = workDir()
        final keep = [d.resolve('S1/Log.final.out').toString(), d.resolve('out.fastq.gz').toString()] as Set
        final cleaner = new WorkDirCleaner(null, false, true)

        when:
        cleaner.clean(d, keep, 'P', 'P (1)', 'ab/cdef')

        then:
        Files.exists(d.resolve('S1/Log.final.out'))
        Files.exists(d.resolve('out.fastq.gz'))
        !Files.exists(d.resolve('S1/_STARtmp'))
        !Files.exists(d.resolve('in.fastq.gz'))

        when:
        cleaner.remove(keep, 'P', 'P (1)', 'ab/cdef')

        then:
        !Files.exists(d.resolve('S1/Log.final.out'))
        !Files.exists(d.resolve('out.fastq.gz'))
    }
}
