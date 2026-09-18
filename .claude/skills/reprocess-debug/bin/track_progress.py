#!/usr/bin/env python3
"""
Project-wide reprocessing progress against the master target list.

"Done" means either of two things, checked independently and unioned:
  - the sample's collection already exists on iRODS under
    <irods-base>/<dataset>/<sample_id> (uploaded by any pipeline, not just
    this one — the master list is not scoped to this repo)
  - this repo's own results/*/starsolo/<dataset>/<sample_id>/ has completed
    STARsolo output (same heuristic as scripts/subset_studies.py:
    output/Gene/ or Log.final.out present)

The join is by sample_id (the target list's first column) alone: for the
iRODS side, against the *leaf* directory name of every two-level-deep
collection under --irods-base, regardless of which dataset folder it sits
under; for the local side, against every results/*/starsolo/*/<id> across
all batches in this repo, not just the batch just debugged. This is a
best-effort join, not a guarantee — see SKILL.md / references/notify.md
for the caveat about non-GEO accessions in the target list that this repo's
pipeline may never touch.

Appends one dated row to --counter-file and prints a short summary (also
available as --json) for the notify step to put in the email/Slack message.
"""
import argparse
import csv
import glob
import json
import os
import subprocess
import sys
import time
from pathlib import Path

DEFAULT_TARGET = "/nfs/cellgeni/projects/reprocessing/irods_datasets/Not_done_hs_10x.sample_table.tsv"
DEFAULT_IRODS_BASE = "/archive/cellgeni/datasets"
DEFAULT_CACHE = "/nfs/cellgeni/reprocessing-runs/.cache/irods-datasets-collections.txt"
DEFAULT_COUNTER = "/nfs/cellgeni/reprocessing-runs/progress-counter.tsv"
CACHE_MAX_AGE_S = 24 * 3600
COUNTER_HEADER = [
    "timestamp", "batch", "run", "total_target", "done_irods", "done_local_only",
    "done_total", "remaining", "pct_done",
]


def read_target_samples(path: str) -> list[str]:
    """First column of the master tab-separated sample table, one id per row."""
    samples = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0].strip():
                samples.append(parts[0].strip())
    return samples


def refresh_irods_cache(cache_path: str, irods_base: str, max_age_s: int) -> bool:
    """
    Regenerate the iRODS collection dump if it is missing or older than
    max_age_s (0 forces a refresh). Requires `module load cellgen/irods` and,
    on this farm, a pre-existing $TMPDIR/<user>/singularity directory for the
    icommands wrapper's Singularity cache — both are the caller's job, not
    this script's, since they are one-off environment setup, not per-run
    state. Returns True if it actually re-ran the query.
    """
    p = Path(cache_path)
    if p.exists() and (max_age_s <= 0 or time.time() - p.stat().st_mtime < max_age_s) and max_age_s > 0:
        return False
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + ".tmp")
    query = f"SELECT COLL_NAME WHERE COLL_NAME like '{irods_base.rstrip('/')}/%'"
    with open(tmp, "w") as out:
        subprocess.run(["iquest", "--no-page", "%s", query], stdout=out, check=True)
    tmp.replace(p)
    return True


def irods_leaf_sample_ids(cache_path: str, irods_base: str) -> set[str]:
    """
    Sample-level leaf names from the cached collection dump: paths exactly
    two segments below irods_base (<base>/<dataset>/<sample>), so a deeper or
    shallower collection (a per-file object listing, or the dataset level
    itself) is not mistaken for a sample.
    """
    prefix = irods_base.rstrip("/") + "/"
    leaves = set()
    with open(cache_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(prefix):
                continue
            parts = line[len(prefix):].split("/")
            if len(parts) == 2 and parts[1]:
                leaves.add(parts[1])
    return leaves


def local_completed_samples(results_glob: str) -> set[str]:
    """Sample dirs under any results/*/starsolo/<dataset>/<sample>/ with STARsolo output."""
    done = set()
    for starsolo_dir in glob.glob(results_glob):
        p = Path(starsolo_dir)
        if not p.is_dir():
            continue
        for dataset_dir in p.iterdir():
            if not dataset_dir.is_dir():
                continue
            for sample_dir in dataset_dir.iterdir():
                if not sample_dir.is_dir():
                    continue
                if (sample_dir / "output" / "Gene").exists() or (sample_dir / "Log.final.out").exists():
                    done.add(sample_dir.name)
    return done


def read_last_counter_row(counter_file: str) -> dict | None:
    if not os.path.exists(counter_file):
        return None
    last = None
    with open(counter_file, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            last = row
    return last


def append_counter_row(counter_file: str, row: dict) -> None:
    exists = os.path.exists(counter_file)
    Path(counter_file).parent.mkdir(parents=True, exist_ok=True)
    with open(counter_file, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COUNTER_HEADER, delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default=DEFAULT_TARGET, help="Master dataset_id/sample_id/.../species table")
    ap.add_argument("--irods-base", default=DEFAULT_IRODS_BASE)
    ap.add_argument("--irods-cache", default=DEFAULT_CACHE)
    ap.add_argument("--cache-max-age", type=int, default=CACHE_MAX_AGE_S,
                     help="Seconds before the iRODS cache is refreshed (default 24h); 0 forces a refresh")
    ap.add_argument("--refresh-irods", action="store_true", help="Shortcut for --cache-max-age 0")
    ap.add_argument("--no-irods", action="store_true",
                     help="Skip the iRODS check entirely (e.g. module not loaded, no cache yet)")
    ap.add_argument("--results-glob", default="results/*/starsolo", help="Relative to cwd (repo root)")
    ap.add_argument("--counter-file", default=DEFAULT_COUNTER)
    ap.add_argument("--batch", default="")
    ap.add_argument("--run", default="")
    ap.add_argument("--json", default=None, help="Also write the summary as JSON to this path")
    ap.add_argument("--dry-run", action="store_true", help="Compute and print, but do not append to --counter-file")
    args = ap.parse_args()

    if args.refresh_irods:
        args.cache_max_age = 0

    targets = read_target_samples(args.target)
    total = len(targets)
    target_set = set(targets)

    irods_note = None
    irods_leaves: set[str] = set()
    if args.no_irods:
        irods_note = "skipped (--no-irods)"
    else:
        try:
            refreshed = refresh_irods_cache(args.irods_cache, args.irods_base, args.cache_max_age)
            irods_leaves = irods_leaf_sample_ids(args.irods_cache, args.irods_base)
            age_h = (time.time() - Path(args.irods_cache).stat().st_mtime) / 3600
            irods_note = f"{'refreshed' if refreshed else 'cached'}, {age_h:.1f}h old, {len(irods_leaves)} sample-level collections"
        except FileNotFoundError:
            irods_note = f"no cache at {args.irods_cache} and could not refresh (iquest not on PATH? `module load cellgen/irods`)"
        except subprocess.CalledProcessError as e:
            irods_note = f"iquest failed (exit {e.returncode}) — is the iRODS session valid? (`iinit`)"

    local_done = local_completed_samples(args.results_glob)

    done_irods = target_set & irods_leaves
    done_local_only = (target_set & local_done) - done_irods
    done_total = done_irods | done_local_only
    remaining = total - len(done_total)
    pct = (100.0 * len(done_total) / total) if total else 0.0

    row = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "batch": args.batch or "",
        "run": args.run or "",
        "total_target": total,
        "done_irods": len(done_irods),
        "done_local_only": len(done_local_only),
        "done_total": len(done_total),
        "remaining": remaining,
        "pct_done": f"{pct:.2f}",
    }

    last = None if args.dry_run else read_last_counter_row(args.counter_file)
    if not args.dry_run:
        append_counter_row(args.counter_file, row)

    delta = None
    if last:
        try:
            delta = row["done_total"] - int(last["done_total"])
        except (KeyError, ValueError):
            delta = None

    summary = {
        **row,
        "target_file": args.target,
        "irods_note": irods_note,
        "delta_done_since_last": delta,
        "previous_measurement": last,
    }

    print(f"target list       : {args.target} ({total} samples)", file=sys.stderr)
    print(f"iRODS              : {irods_note}", file=sys.stderr)
    print(f"done on iRODS      : {len(done_irods)}", file=sys.stderr)
    print(f"done locally only  : {len(done_local_only)}  (STARsolo complete, not yet seen on iRODS)", file=sys.stderr)
    print(f"done total         : {len(done_total)} / {total}  ({pct:.2f}%)", file=sys.stderr)
    print(f"remaining          : {remaining}", file=sys.stderr)
    if delta is not None:
        sign = "+" if delta >= 0 else ""
        print(f"change since last measurement ({last['timestamp']}): {sign}{delta}", file=sys.stderr)
    if not args.dry_run:
        print(f"appended to {args.counter_file}", file=sys.stderr)

    if args.json:
        with open(args.json, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"wrote {args.json}", file=sys.stderr)

    # also on stdout, for a caller that isn't reading stderr
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
