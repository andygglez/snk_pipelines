#!/usr/bin/env python3
"""Command-line front end for the Snakemake ports in this repository.

    pipeline list                              pipelines and their image downloads
    pipeline config <name> [-p PATH]           copy the default config.yaml
    pipeline pull   <name> [-c CONFIG]         download missing Singularity images
    pipeline run    <name> -c CONFIG [-t N]    run the port with apptainer

<name> is the pipeline's official name (rnaseq, chipseq, atacseq, hic, wgbs);
the port directory's suffix also works (methylseq -> wgbs). Runs happen in
the current directory (or --workdir): outputs and relative paths in CONFIG are
resolved there, not inside the repository.

Images: every rule runs from the local image listed under `containers:` in
snk-<name>/config/containers.yaml. Snakemake never downloads them, and compute
nodes are offline, so `pipeline pull` them on a machine with internet (e.g. a
Rorqual login node) first. `pipeline run` refuses to start while any image is
missing.

Install: ln -s <repo>/pipeline.py ~/.local/bin/pipeline
"""

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import urllib.request

try:
    import yaml
except ImportError:
    sys.exit("pipeline: PyYAML is required (pip install pyyaml, or use the "
             "Python of the environment that has Snakemake).")

REPO = os.path.dirname(os.path.realpath(__file__))
PREFIX = "snk-"
# Official pipeline name -> port directory, in `pipeline list` order. Ports not
# listed here are shown after these under their directory suffix.
NAMES = {
    "rnaseq": "snk-rnaseq",
    "chipseq": "snk-chipseq",
    "atacseq": "snk-atacseq",
    "hic": "snk-hic",
    "wgbs": "snk-methylseq",
}


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def die(msg):
    sys.exit(f"pipeline: {msg}")


def port_dirs():
    """Official name -> port directory, for every port that has a Snakefile."""
    found = {d for d in os.listdir(REPO) if d.startswith(PREFIX)
             and os.path.isfile(os.path.join(REPO, d, "Snakefile"))}
    out = {name: d for name, d in NAMES.items() if d in found}
    for d in sorted(found - set(NAMES.values())):
        out[d[len(PREFIX):]] = d
    return out


def ports():
    return list(port_dirs())


def port_dir(name):
    dirs = port_dirs()
    d = dirs.get(name)
    if d is None and PREFIX + name in dirs.values():  # e.g. methylseq
        d = PREFIX + name
    if d is None:
        die(f"unknown pipeline {name!r}; available: {', '.join(dirs)}")
    return os.path.join(REPO, d)


def load_yaml(path):
    try:
        with open(path) as fh:
            return yaml.safe_load(fh) or {}
    except FileNotFoundError:
        die(f"file not found: {path}")


def images(name, user_config=None):
    """Rule key -> (absolute image path, source URL), with the user config's
    `containers:` overrides applied the way Snakemake merges config files."""
    pdir = port_dir(name)
    data = load_yaml(os.path.join(pdir, "config", "containers.yaml"))
    paths = dict(data.get("containers", {}))
    if user_config:
        paths.update(load_yaml(user_config).get("containers") or {})
    urls = data.get("container_urls", {})
    out = {}
    for key, path in paths.items():
        if not os.path.isabs(path):
            path = os.path.join(pdir, path)
        out[key] = (path, urls.get(key))
    return out


def missing_images(name, user_config=None):
    return {k: v for k, v in images(name, user_config).items()
            if not os.path.isfile(v[0])}


def auto_binds(configfile, workdir):
    """Directories the containers need to see the paths in the user config.

    Every string value of the config (comma-separated lists too) that names an
    existing path, or an absolute path whose parent exists (outdir, workdir), is
    bound by its directory. Small CSV/TSV/TXT files (samplesheets, BBSplit and
    rRNA lists) are scanned for paths the same way. Both the path as written
    and its symlink-resolved target are bound: apptainer starts with --home set
    to the run directory, so neither $HOME nor /project symlink targets exist
    inside the container otherwise."""
    dirs, seen = [], set()

    def add_dir(d):
        for p in (os.path.abspath(d), os.path.realpath(d)):
            if p != os.sep and "," not in p and ":" not in p:
                dirs.append(p)

    def visit(token, base):
        token = token.strip().strip('"\'')
        if not token or "://" in token or "\n" in token or len(token) > 4096:
            return
        path = os.path.abspath(os.path.join(base, os.path.expanduser(token)))
        if path in seen:
            return
        seen.add(path)
        if os.path.isdir(path):
            add_dir(path)
        elif os.path.exists(path):
            add_dir(os.path.dirname(path))
            real = os.path.realpath(path)
            if real != path:
                add_dir(os.path.dirname(real))
            if (path.endswith((".csv", ".tsv", ".txt"))
                    and os.path.getsize(path) < (1 << 20)):
                with open(path, errors="replace") as fh:
                    for line in fh:
                        for cell in line.replace("\t", ",").split(","):
                            visit(cell, os.path.dirname(path))
        elif os.path.isabs(token) and os.path.isdir(os.path.dirname(path)):
            add_dir(os.path.dirname(path))  # output dir not created yet

    def walk(value):
        if isinstance(value, dict):
            for k, v in value.items():
                if k not in ("containers", "container_urls"):
                    walk(v)
        elif isinstance(value, list):
            for v in value:
                walk(v)
        elif isinstance(value, str):
            for token in value.split(","):
                visit(token, workdir)

    walk(load_yaml(configfile))
    # Drop directories already covered by a bound parent.
    out = []
    for d in sorted(set(dirs)):
        if not any(d.startswith(p.rstrip(os.sep) + os.sep) for p in out):
            out.append(d)
    return out


def find_snakemake(override):
    # The snakemake installed next to this Python (uv sync / uv tool install)
    # comes before PATH: `uv tool` exposes only `pipeline`, not snakemake.
    sibling = os.path.join(os.path.dirname(sys.executable), "snakemake")
    exe = (override or os.environ.get("PIPELINE_SNAKEMAKE")
           or (sibling if os.access(sibling, os.X_OK) else None)
           or shutil.which("snakemake"))
    if not exe:
        die("snakemake not found on PATH. Activate its environment, or pass "
            "--snakemake / set PIPELINE_SNAKEMAKE to the executable.")
    return exe


# -----------------------------------------------------------------------------
# Subcommands
# -----------------------------------------------------------------------------
def cmd_list(args):
    color = sys.stdout.isatty() and not os.environ.get("NO_COLOR")

    def paint(text, code):
        return f"\033[{code}m{text}\033[0m" if color else text

    width = 32
    for i, name in enumerate(ports()):
        total = len(images(name))
        have = total - len(missing_images(name))
        filled = round(width * have / total) if total else 0
        tint = "32" if have == total else "33" if have else "31"  # green/yellow/red
        bar = paint("█" * filled, tint) + paint("░" * (width - filled), "2")
        pct = 100 * have // total if total else 0
        if i:
            print()
        print(paint(name, "1"))
        print(f"  {bar}  {have:>2}/{total:<2} images  {paint(f'{pct:>3}%', tint)}")


def cmd_config(args):
    src = os.path.join(port_dir(args.name), "config", "config.yaml")
    dest = args.path
    if os.path.isdir(dest) or dest.endswith(os.sep):
        dest = os.path.join(dest, f"{args.name}.config.yaml")
    if os.path.exists(dest) and not args.force:
        die(f"{dest} exists; pass --force to overwrite")
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    shutil.copyfile(src, dest)
    print(f"Wrote {dest}")
    print("Edit at least input, outdir and the reference paths. Relative paths "
          "are resolved from the directory you run `pipeline run` in.")


def cmd_pull(args):
    miss = missing_images(args.name, args.configfile)
    if not miss:
        print(f"{args.name}: all images present")
        return
    failed = []
    # One at a time, on purpose: parallel pulls overload small machines.
    for i, (key, (path, url)) in enumerate(miss.items(), 1):
        if not url:
            print(f"[{i}/{len(miss)}] {key}: no URL in container_urls; skipped")
            failed.append(key)
            continue
        print(f"[{i}/{len(miss)}] {key}: {url}", flush=True)
        if args.dry_run:
            continue
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if url.startswith("docker://"):
            # Docker-only images (e.g. Parabricks) are converted by apptainer.
            rc = subprocess.call(["apptainer", "pull", path, url])
            if rc != 0:
                print(f"    failed: apptainer pull exited with {rc}")
                failed.append(key)
            continue
        part = path + ".part"
        try:
            # community-cr-prod.seqera.io (Cloudflare) answers 403 to urllib's
            # default "Python-urllib/3.x" User-Agent, so send our own.
            req = urllib.request.Request(url, headers={"User-Agent": "pipeline/1.0"})
            with urllib.request.urlopen(req, timeout=60) as resp, open(part, "wb") as fh:
                shutil.copyfileobj(resp, fh, length=1 << 20)
            os.replace(part, path)
            print(f"    -> {path} ({os.path.getsize(path) / 1e6:.0f} MB)")
        except Exception as exc:  # network errors, HTTP errors, disk full
            if os.path.exists(part):
                os.remove(part)
            print(f"    failed: {exc}")
            failed.append(key)
    if failed:
        die(f"{len(failed)} image(s) not downloaded: {', '.join(failed)}")


def cmd_run(args):
    pdir = port_dir(args.name)
    configfile = os.path.abspath(args.configfile)
    if not os.path.isfile(configfile):
        die(f"config file not found: {args.configfile}")
    workdir = os.path.abspath(args.workdir)

    miss = missing_images(args.name, configfile)
    if miss:
        msg = (f"{len(miss)} image(s) missing: {', '.join(miss)}. Run "
               f"`pipeline pull {args.name}` on a machine with internet first.")
        if args.dry_run or args.no_check:
            print(f"pipeline: warning: {msg}", file=sys.stderr)
        else:
            die(msg)

    binds = [workdir] + [os.path.abspath(b) for b in args.bind]
    if not args.no_auto_bind:
        binds += auto_binds(configfile, workdir)
    cmd = [
        find_snakemake(args.snakemake),
        "--snakefile", os.path.join(pdir, "Snakefile"),
        "--directory", workdir,
        # Port defaults first, the user's config last (later files win). Passed
        # explicitly because the Snakefile's own `configfile:` paths are
        # relative and would be looked up in workdir.
        "--configfile",
        os.path.join(pdir, "config", "config.yaml"),
        os.path.join(pdir, "config", "containers.yaml"),
        configfile,
        "--cores", str(args.threads),
        "--sdm", "apptainer",
        "--apptainer-args", "--bind " + ",".join(dict.fromkeys(binds)),
    ]
    if args.dry_run:
        cmd.append("--dry-run")
    cmd += args.snakemake_args

    if not args.dry_run and not (shutil.which("apptainer") or shutil.which("singularity")):
        print("pipeline: warning: apptainer not on PATH (on Rorqual: module load "
              "apptainer)", file=sys.stderr)
    print("+ " + shlex.join(cmd), file=sys.stderr, flush=True)
    sys.exit(subprocess.call(cmd))


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        prog="pipeline", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("list", help="list the pipelines and their image downloads")
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("config", help="copy a pipeline's default config.yaml")
    p.add_argument("name")
    p.add_argument("-p", "--path", default=".",
                   help="directory (writes <name>.config.yaml) or file path [.]")
    p.add_argument("-f", "--force", action="store_true", help="overwrite")
    p.set_defaults(func=cmd_config)

    p = sub.add_parser("pull", help="download missing images")
    p.add_argument("name")
    p.add_argument("-c", "--configfile",
                   help="user config, if it overrides `containers:`")
    p.add_argument("-n", "--dry-run", action="store_true",
                   help="only list what would be downloaded")
    p.set_defaults(func=cmd_pull)

    p = sub.add_parser(
        "run", help="run a pipeline",
        description="Extra Snakemake options go after `--`, e.g. "
                    "pipeline run rnaseq -c my.yaml -- --rerun-incomplete")
    p.add_argument("name")
    p.add_argument("-c", "--configfile", required=True, help="your config YAML")
    p.add_argument("-t", "--threads", type=int, default=1, help="cores [1]")
    p.add_argument("-d", "--workdir", default=".",
                   help="directory to run in (outputs, relative paths) [.]")
    p.add_argument("-B", "--bind", action="append", default=[],
                   help="extra path to bind into the containers (repeatable); "
                        "the workdir and the directories of every path in the "
                        "config and its samplesheets are bound automatically")
    p.add_argument("--no-auto-bind", action="store_true",
                   help="bind only the workdir and -B paths")
    p.add_argument("-n", "--dry-run", action="store_true")
    p.add_argument("--no-check", action="store_true",
                   help="start even if images are missing")
    p.add_argument("--snakemake", help="snakemake executable [PATH lookup]")
    p.set_defaults(func=cmd_run)

    argv = sys.argv[1:]
    extra = []
    if "--" in argv:
        i = argv.index("--")
        argv, extra = argv[:i], argv[i + 1:]
    args = parser.parse_args(argv)
    args.snakemake_args = extra
    args.func(args)


if __name__ == "__main__":
    main()
