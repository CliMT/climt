# scripts/build_experiments.py
"""Regenerate experiment/chapter artifacts declared in docs/**/sources.yml.

Each sources.yml lists artifacts; each artifact maps an output figure (and
optional captured-stdout file) to the command that produces it and the source
files whose content determines its validity. Deps are content-hashed (sha256,
glob-expanded) and compared to a stored _artifacts/.hashes.json; only stale
artifacts (output missing or dep-hash changed) are re-run. Hashes, not mtimes,
so branch switches and network filesystems don't trigger spurious rebuilds.

    python scripts/build_experiments.py            # regenerate stale artifacts
    python scripts/build_experiments.py --check    # CI: report + exit 1 if stale

``--only`` narrows the run to artifacts whose repo-relative output path matches a
glob, which matters because a full sweep re-runs multi-hundred-day RCE
integrations and takes over an hour. A change confined to one component usually
invalidates one manifest's worth of artifacts, so regenerate just those:

    python scripts/build_experiments.py --only 'docs/experiments/*cork-co2-bands/**'
    python scripts/build_experiments.py --only '**/throughput.npz' --only '**/throughput.png'

``--only`` selects; it does not follow dependencies. An artifact that consumes
another artifact's output (a figure reading an .npz) is only rebuilt if it is
itself selected, so include the consumers when you narrow to a producer.

An artifact carrying ``manual: true`` is skipped everywhere — never reported by
``--check``, never rebuilt — unless ``--manual`` is passed. That is for outputs
whose content depends on the machine rather than on the code, wall-clock
benchmarks above all: hash-gating a timing means any edit to a dep forces the
numbers to be re-recorded on whoever's laptop runs the sweep, so hardware noise
lands in the repo as a real diff and a stale-artifact CI failure. Regenerate one
deliberately, on a machine whose numbers you would stand behind:

    python scripts/build_experiments.py --manual --only '**/throughput.*'
"""
import argparse
import fnmatch
import glob
import hashlib
import json
import os
import subprocess
import sys

import yaml


def _expand_deps(repo_root, deps):
    paths = []
    for pattern in deps:
        matches = sorted(glob.glob(os.path.join(repo_root, pattern), recursive=True))
        if not matches:
            raise FileNotFoundError(f"dep pattern matched no files: {pattern}")
        paths.extend(matches)
    return paths


def _hash_files(repo_root, paths):
    # Hash the repo-root-relative path (not the absolute path) plus file bytes,
    # so the digest is stable across machines and git worktrees: the same file
    # tree yields the same hash regardless of where it is checked out.
    h = hashlib.sha256()
    for p in sorted(set(paths)):
        h.update(os.path.relpath(p, repo_root).encode())
        with open(p, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 16), b""):
                h.update(chunk)
    return h.hexdigest()


def _load_hashes(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def _stale(repo_root, manifest_dir, artifact, stored):
    out = os.path.join(manifest_dir, artifact["out"])
    current = _hash_files(repo_root, _expand_deps(repo_root, artifact["deps"]))
    if not os.path.exists(out):
        return True, current
    return (stored.get(artifact["out"]) != current), current


def _selected(repo_root, manifest_dir, artifact, patterns):
    """Match --only globs against the artifact's repo-relative output path."""
    if not patterns:
        return True
    rel = os.path.relpath(
        os.path.join(manifest_dir, artifact["out"]), repo_root)
    return any(fnmatch.fnmatch(rel, p) for p in patterns)


def _regenerate_manifest(repo_root, manifest_path, check_only, only, manual):
    manifest_dir = os.path.dirname(manifest_path)
    with open(manifest_path) as f:
        manifest = yaml.safe_load(f)
    hashes_path = os.path.join(manifest_dir, "_artifacts", ".hashes.json")
    os.makedirs(os.path.dirname(hashes_path), exist_ok=True)
    stored = _load_hashes(hashes_path)
    any_stale = False
    n_selected = 0
    # Artifacts run in declaration order: producers must precede consumers within
    # a manifest (a figure that reads another artifact's output lists it later).
    for art in manifest["artifacts"]:
        if not _selected(repo_root, manifest_dir, art, only):
            continue
        n_selected += 1
        if art.get("manual", False) and not manual:
            # Counted as selected above, so an --only naming just this artifact
            # reports why nothing happened instead of "matched no artifacts".
            if only:
                print(f"skipping {art['out']} (manual; pass --manual to run it)")
            continue
        stale, current = _stale(repo_root, manifest_dir, art, stored)
        if not stale:
            continue
        any_stale = True
        if check_only:
            print(f"STALE: {art['out']} (dep changed or output missing)")
            continue
        print(f"regenerating {art['out']} ...")
        done = subprocess.run(art["cmd"], shell=True, cwd=repo_root,
                              capture_output=True, text=True)
        if done.returncode != 0:
            # Emit both streams to stderr so CI captures the full failure context.
            sys.stderr.write(done.stdout + done.stderr)
            raise RuntimeError(f"cmd failed for {art['out']}: {art['cmd']}")
        if "out_txt" in art:
            with open(os.path.join(manifest_dir, art["out_txt"]), "w") as f:
                f.write(done.stdout)
        if not os.path.exists(os.path.join(manifest_dir, art["out"])):
            raise RuntimeError(
                f"cmd for {art['out']} did not write the output file")
        stored[art["out"]] = current
    # `stored` is only ever added to, so a narrowed run preserves the recorded
    # hashes of the artifacts it skipped.
    if not check_only and n_selected:
        with open(hashes_path, "w") as f:
            json.dump(stored, f, indent=2, sort_keys=True)
    return any_stale, n_selected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="report staleness and exit 1 if any; no execution")
    ap.add_argument("--root", default=os.getcwd(), help="repo root (default cwd)")
    ap.add_argument("--only", action="append", metavar="GLOB",
                    help="restrict to artifacts whose repo-relative output path "
                         "matches GLOB; repeatable. Selects only — consumers of a "
                         "selected artifact are not pulled in automatically.")
    ap.add_argument("--manual", action="store_true",
                    help="also process artifacts marked `manual: true` in their "
                         "manifest, which are otherwise never checked or rebuilt.")
    args = ap.parse_args()
    manifests = sorted(glob.glob(
        os.path.join(args.root, "docs", "**", "sources.yml"), recursive=True))
    any_stale = False
    total_selected = 0
    for m in manifests:
        stale, n = _regenerate_manifest(
            args.root, m, args.check, args.only, args.manual)
        any_stale |= stale
        total_selected += n
    if args.only and not total_selected:
        sys.stderr.write(
            f"--only matched no artifacts: {args.only}\n"
            "Patterns match the repo-relative output path, e.g. "
            "'docs/experiments/*cork-co2-bands/**'.\n")
        sys.exit(2)
    if args.check and any_stale:
        sys.stderr.write("Stale artifacts. Run `make experiments` and commit.\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
