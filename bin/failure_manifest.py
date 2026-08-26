#!/usr/bin/env python3
"""Write a failure manifest for a Nextflow run.

The pipeline sets errorStrategy 'ignore' almost everywhere, so failures never
reach the exit status and nothing records them. Post-mortems have had to be mined
out of the driver's console log by hand. This collects the same information
directly from the work dirs, which persist.

It is driven from `workflow.onComplete` in nextflow.config, and is also safe to
run by hand against any past run:

    bin/failure_manifest.py --workdir nf-work --out reports/failures.tsv

Work dirs are the source of truth rather than the trace file: every task attempt
leaves a `.exitcode` and a `.command.run` carrying its process name and tag, so
the manifest does not depend on when Nextflow flushes its trace. Pass --since to
restrict the scan to one run, since nf-work accumulates across runs.

A note on counting, because two units are easy to conflate. Each row is one
failed task *instance*. A single run accession can own several instances -- one
per URL, and one per dataset holding its sample -- and in Nextflow those tasks
can share an identical tag. Counting by tag collapses them: in batch 2 that
undercounted permanent failures by 8. The summary reports both units.
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict

NAME_RE = re.compile(r"^### name: '(?P<name>.*)'\s*$")

# Lines worth surfacing as "the error". A real message almost always announces
# itself with one of these.
ERROR_RE = re.compile(
    r"^\s*(?:\[ERROR\s*\]|\[ERROR\]|ERROR:|ERROR\b|FATAL|Traceback \(most recent call last\)|"
    r"[A-Za-z_.]*(?:Error|Exception):|usage:|EXITING because of)"
)
# Network trouble is the real cause for a download task, but it repeats hundreds of
# times per log, so it is a fallback rather than a first choice.
NETWORK_RE = re.compile(r"Read error at byte|Giving up\.|unable to resolve host|"
                        r"Connection timed out|failed: Connection|No route to host")
# wget's progress bar. Never informative, and there can be hundreds of thousands.
PROGRESS_RE = re.compile(r"^\s*\d+K[\s.]|^\s*\d+%\s|\.\.\.\.\.\.")


def task_name(workdir):
    """Return the Nextflow task name recorded in .command.run, or None."""
    path = os.path.join(workdir, ".command.run")
    try:
        with open(path, errors="replace") as fh:
            for i, line in enumerate(fh):
                if i > 40:
                    break
                m = NAME_RE.match(line)
                if m:
                    return m.group("name")
    except OSError:
        pass
    return None


def first_error(workdir, limit=400):
    """Best-effort one-line summary of why a task failed.

    Preference order: an explicit error message, then a network failure, then any
    substantial line. The tiers matter because a download that times out has no
    ERROR: line at all -- its cause is a "Read error ... Giving up." buried in
    progress output -- and reporting nothing for those made 11 of batch 2's rows
    blank, which is precisely the case that most needed a message.
    """
    best = {"error": None, "network": None, "other": None}
    for fname in (".command.err", ".command.log"):
        path = os.path.join(workdir, fname)
        try:
            with open(path, errors="replace") as fh:
                for line in fh:
                    line = line.rstrip("\n").strip()
                    if not line or PROGRESS_RE.search(line):
                        continue
                    if ERROR_RE.match(line):
                        if best["error"] is None:
                            best["error"] = line[:limit]
                        # An explicit error is the best we will do; stop early.
                        return best["error"]
                    if NETWORK_RE.search(line):
                        if best["network"] is None:
                            best["network"] = line[:limit]
                    elif best["other"] is None and len(line) > 8:
                        best["other"] = line[:limit]
        except OSError:
            continue
    return best["error"] or best["network"] or best["other"] or ""


def exit_code(workdir):
    try:
        with open(os.path.join(workdir, ".exitcode")) as fh:
            return fh.read().strip()
    except OSError:
        return None


def is_array_driver(workdir):
    """True for an LSF array-job dispatcher rather than a real task.

    Processes with `array = N` set (WGET10X uses 20) get an extra work dir per
    batch whose .command.sh only execs the N child work dirs. It carries a numeric
    tag and a non-zero .exitcode whenever a child fails, so counting it reports
    failures that do not exist -- 11 phantom WGET10X rows on the batch 2 data,
    with no error text because the real output goes to the children's logs. The
    children are scanned in their own right, so skip the dispatcher.
    """
    try:
        with open(os.path.join(workdir, ".command.sh"), errors="replace") as fh:
            return "nxf_array_task_dir" in fh.read(4096)
    except OSError:
        return False


def scan(workroot, since=None):
    """Yield (workdir, name, exit) for every task attempt found."""
    if not os.path.isdir(workroot):
        return
    for bucket in sorted(os.listdir(workroot)):
        bpath = os.path.join(workroot, bucket)
        if len(bucket) != 2 or not os.path.isdir(bpath):
            continue
        try:
            entries = sorted(os.listdir(bpath))
        except OSError:
            continue
        for h in entries:
            wd = os.path.join(bpath, h)
            code = exit_code(wd)
            if code is None:
                continue
            if since is not None:
                try:
                    if os.path.getmtime(os.path.join(wd, ".exitcode")) < since:
                        continue
                except OSError:
                    continue
            if is_array_driver(wd):
                continue
            yield wd, task_name(wd) or "UNKNOWN", code


def split_name(name):
    """'A:B:PROC (tag)' -> ('PROC', 'tag')"""
    m = re.match(r"^(?P<proc>[^ ]+?)(?: \((?P<tag>.*)\))?$", name)
    if not m:
        return name, ""
    proc = m.group("proc").split(":")[-1]
    return proc, m.group("tag") or ""


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--workdir", default="nf-work", help="Nextflow work directory")
    p.add_argument("--out", default="reports/failures.tsv", help="manifest to write")
    p.add_argument("--since", type=float, default=None,
                   help="only tasks whose .exitcode is newer than this epoch time")
    p.add_argument("--quiet", action="store_true", help="write the file, print nothing")
    args = p.parse_args()

    ok_names = set()
    failed = defaultdict(list)          # name -> [(mtime, workdir, exit), ...]
    for wd, name, code in scan(args.workdir, args.since):
        if code == "0":
            ok_names.add(name)
        else:
            try:
                mtime = os.path.getmtime(os.path.join(wd, ".exitcode"))
            except OSError:
                mtime = 0.0
            failed[name].append((mtime, wd, code))

    # One row per task, not per attempt. Every retry leaves its own work dir with a
    # non-zero .exitcode, so counting those counts attempts: batch 2's 208 permanently
    # failed RENAME10XRUN tasks leave 416 of them. The last attempt is the one whose
    # error matters, so report that and carry the attempt count alongside.
    rows = []
    for name, attempts in failed.items():
        attempts.sort()
        _, last_wd, last_code = attempts[-1]
        proc, tag = split_name(name)
        rows.append({
            "process": proc,
            "tag": tag,
            "exit": last_code,
            "attempts_failed": len(attempts),
            # A task whose name also appears with exitcode 0 was retried and
            # recovered. "Usually": concurrent tasks can share a tag, so for those
            # this is a hint rather than a verdict. The tag now carries the file
            # name for WGET10X, which removes the common case of that ambiguity.
            "recovered": "yes" if name in ok_names else "no",
            "workdir": os.path.abspath(last_wd),
            "first_error": first_error(last_wd),
        })

    rows.sort(key=lambda r: (r["recovered"], r["process"], r["tag"]))
    cols = ["process", "tag", "exit", "attempts_failed", "recovered", "workdir", "first_error"]

    outdir = os.path.dirname(os.path.abspath(args.out))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

    if args.quiet:
        return 0

    hard = [r for r in rows if r["recovered"] == "no"]
    soft = [r for r in rows if r["recovered"] == "yes"]
    attempts = sum(r["attempts_failed"] for r in rows)
    print(f"Failure manifest: {args.out}")
    print(f"  tasks that failed     : {len(rows)}  ({attempts} failed attempts incl. retries)")
    print(f"  never succeeded       : {len(hard)}   <- what the run actually lost")
    print(f"  recovered on retry    : {len(soft)}")
    if hard:
        print("\n  never succeeded, by process and exit code:")
        for (proc, code), n in Counter((r["process"], r["exit"]) for r in hard).most_common():
            print(f"    {proc:<28} exit {code:<5} {n}")
        print("\n  most common first error line:")
        norm = defaultdict(int)
        for r in hard:
            key = re.sub(r"[SED]RR\d+|GSM\d+|GSE\d+", "<ACC>", r["first_error"])[:110]
            norm[key] += 1
        for msg, n in sorted(norm.items(), key=lambda x: -x[1])[:8]:
            print(f"    {n:>5}  {msg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
