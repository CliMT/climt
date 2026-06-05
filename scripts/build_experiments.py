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
"""
import argparse
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


def _regenerate_manifest(repo_root, manifest_path, check_only):
    manifest_dir = os.path.dirname(manifest_path)
    with open(manifest_path) as f:
        manifest = yaml.safe_load(f)
    hashes_path = os.path.join(manifest_dir, "_artifacts", ".hashes.json")
    os.makedirs(os.path.dirname(hashes_path), exist_ok=True)
    stored = _load_hashes(hashes_path)
    any_stale = False
    # Artifacts run in declaration order: producers must precede consumers within
    # a manifest (a figure that reads another artifact's output lists it later).
    for art in manifest["artifacts"]:
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
    if not check_only:
        with open(hashes_path, "w") as f:
            json.dump(stored, f, indent=2, sort_keys=True)
    return any_stale


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="report staleness and exit 1 if any; no execution")
    ap.add_argument("--root", default=os.getcwd(), help="repo root (default cwd)")
    args = ap.parse_args()
    manifests = sorted(glob.glob(
        os.path.join(args.root, "docs", "**", "sources.yml"), recursive=True))
    any_stale = False
    for m in manifests:
        any_stale |= _regenerate_manifest(args.root, m, args.check)
    if args.check and any_stale:
        sys.stderr.write("Stale artifacts. Run `make experiments` and commit.\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
