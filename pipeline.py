#!/usr/bin/env python3
"""Command-line front end for the Snakemake ports in this repository.

    pipeline list                              pipelines and their image downloads
    pipeline config <name> [-p PATH]           copy the default config.yaml
    pipeline pull   <name> [-c CONFIG]         download missing Singularity images
    pipeline run    <name> -c CONFIG [-t N]    run the port with apptainer
    pipeline references -p DIR [-a ASM] [--pipeline NAME ...]
                                               download/build reference files

<name> is the pipeline's official name (rnaseq, chipseq, atacseq, hic, wgbs);
the port directory's suffix also works (methylseq -> wgbs). Runs happen in
the current directory (or --workdir): outputs and relative paths in CONFIG are
resolved there, not inside the repository.

Images: every rule runs from the local image listed under `containers:` in
snk-<name>/config/containers.yaml. Snakemake never downloads them, and compute
nodes are offline, so `pipeline pull` them on a machine with internet (e.g. a
Rorqual login node) first. `pipeline run` refuses to start while any image is
missing.

References: `pipeline references` fills DIR/<assembly>/ with the base files of
one assembly from references/assemblies.yaml (FASTA, GTF, blacklist, cytoBand)
and the light files each port derives from them (snk-<name>/config/
references.yaml), built in the port's images. It builds no aligner index and
the pipelines never read DIR by themselves: copy the paths you need from
DIR/<assembly>/README.tsv into your run config. Needs internet and apptainer
(a Rorqual login node), and the images (`pipeline pull <name>`).

Install: ln -s <repo>/pipeline.py ~/.local/bin/pipeline
"""

import argparse
import gzip
import itertools
import os
import re
import shlex
import shutil
import string
import subprocess
import sys
import urllib.parse
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
    """Rule key -> (absolute image path, source), with the user config's
    `containers:` overrides applied the way Snakemake merges config files.

    The source is the URL from `container_urls:`, or "def:<path>" for an image
    built from an Apptainer definition file (`container_defs:`)."""
    pdir = port_dir(name)
    data = load_yaml(os.path.join(pdir, "config", "containers.yaml"))
    paths = dict(data.get("containers", {}))
    if user_config:
        paths.update(load_yaml(user_config).get("containers") or {})
    urls = dict(data.get("container_urls", {}))
    for key, deffile in (data.get("container_defs") or {}).items():
        urls[key] = "def:" + (deffile if os.path.isabs(deffile) else os.path.join(pdir, deffile))
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
                if k not in ("containers", "container_urls", "container_defs"):
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
            print(f"[{i}/{len(miss)}] {key}: no URL in container_urls or container_defs; skipped")
            failed.append(key)
            continue
        if os.path.isfile(path):  # shared with a key pulled earlier in this loop
            print(f"[{i}/{len(miss)}] {key}: already pulled ({os.path.basename(path)})")
            continue
        print(f"[{i}/{len(miss)}] {key}: {url}", flush=True)
        if args.dry_run:
            continue
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if url.startswith("def:"):
            # Port-specific image (no nf-core module): built from its definition file.
            rc = subprocess.call(["apptainer", "build", "--fakeroot", path, url[len("def:"):]])
            if rc != 0:
                print(f"    failed: apptainer build exited with {rc}")
                failed.append(key)
            continue
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


def official_name(name):
    """The `pipeline list` name of a port given by name or directory suffix."""
    pdir = os.path.basename(port_dir(name))
    return next(n for n, d in port_dirs().items() if d == pdir)


def choose_assembly(catalog):
    names = list(catalog)
    if not sys.stdin.isatty():
        die(f"--assembly is required without a terminal; available: {', '.join(names)}")
    print("Available assemblies:")
    for i, name in enumerate(names, 1):
        print(f"  {i}) {name:6s} {catalog[name].get('description', '')}")
    while True:
        try:
            ans = input(f"Assembly [1-{len(names)} or name]: ").strip()
        except EOFError:
            die("no assembly chosen")
        if ans.isdigit() and 1 <= int(ans) <= len(names):
            return names[int(ans) - 1]
        if ans in catalog:
            return ans
        print(f"  not an option: {ans!r}")


def apptainer_runner(image, cmd, cwd, binds):
    """Run `cmd` with bash inside `image`, working directory `cwd`."""
    exe = shutil.which("apptainer") or shutil.which("singularity")
    binds = list(dict.fromkeys(binds + [os.path.realpath(b) for b in binds]))
    env = dict(os.environ, APPTAINERENV_MPLCONFIGDIR="/tmp/mplcfg",
               SINGULARITYENV_MPLCONFIGDIR="/tmp/mplcfg")
    # Lmod exports a `which` shell function that breaks inside containers.
    cmd = "unset -f which 2>/dev/null || true\n" + cmd
    return subprocess.call([exe, "exec", "--bind", ",".join(binds), "--pwd", cwd,
                            image, "bash", "-euo", "pipefail", "-c", cmd], env=env)


def cmd_references(args):
    try:
        catalog = load_catalog(os.path.join(REPO, "references", "assemblies.yaml"))
        assembly = args.assembly or choose_assembly(catalog)
        if assembly not in catalog:
            die(f"unknown assembly {assembly!r}; available: {', '.join(catalog)}")
        names = [official_name(n) for n in args.pipeline] if args.pipeline else ports()
        loaded = []
        for name in dict.fromkeys(names):
            pdir = port_dir(name)
            if not has_spec(pdir):
                if args.pipeline:
                    die(f"{name}: no config/references.yaml in {pdir}")
                continue
            imgs = {k: v[0] for k, v in images(name).items()}
            loaded.append(load_port(name, pdir, imgs))
        plan = plan_refs(catalog, assembly, loaded, REPO)
    except RefError as exc:
        die(str(exc))

    adir = os.path.join(os.path.abspath(os.path.expanduser(args.path)), assembly)
    pending = plan.pending(adir, args.force)
    builds = [j for j in pending if j.kind == "build"]
    missing = sorted({(j.port["name"], j.container) for j in builds if not os.path.isfile(j.image)})
    print(f"{assembly} -> {adir} ({', '.join(p['name'] for p in loaded)})")

    if args.dry_run:
        todo = {id(j) for j in pending}
        for job in plan.jobs:
            if job.unavailable:
                what = f"unavailable: {job.unavailable}"
            elif id(job) not in todo:
                what = "exists"
            elif job.kind == "fetch":
                what = f"download {job.source}"
            else:
                what = f"build in {job.container}" + (" (image missing)" if not os.path.isfile(job.image) else "")
            print(f"  {job.name}: {', '.join(job.outputs)}  [{what}]")
        if missing:
            print(f"pipeline: warning: {len(missing)} image(s) missing; run "
                  + "; ".join(sorted({f'pipeline pull {p}' for p, _ in missing})), file=sys.stderr)
        return

    if builds and not (shutil.which("apptainer") or shutil.which("singularity")):
        die("apptainer is not on PATH (on Rorqual: module load apptainer)")
    if missing:
        die(f"{len(missing)} image(s) missing ({', '.join(f'{p}:{k}' for p, k in missing)}); run "
            + "; ".join(sorted({f"`pipeline pull {p}`" for p, _ in missing})) + " first")

    status = execute_refs(plan, os.path.dirname(adir), apptainer_runner, force=args.force)
    print(f"{assembly}: {refs_summary(status)}. File list: {os.path.join(adir, 'README.tsv')}")
    if any(s in ("failed", "blocked") for s in status.values()):
        sys.exit(1)


# -----------------------------------------------------------------------------
# References engine (`pipeline references`)
#
# The shared catalog (references/assemblies.yaml) says where each base file of an
# assembly comes from. Each port's config/references.yaml lists the files that
# port needs: a *fetch* entry copies a catalog file, and a *build* entry runs a
# shell command, in one of the port's images, on files made by other entries.
# Everything lands flat in <refdir>/<assembly>/ (see the design spec,
# docs/superpowers/specs/2026-11-06-pipeline-references-design.md). The engine
# takes a runner that executes a command in an image (apptainer_runner, or
# local_runner in the tests), so it can be tested without containers.
# -----------------------------------------------------------------------------
FETCH_KEYS = {"catalog", "out"}
BUILD_KEYS = {"from", "out", "container", "shell", "params", "foreach", "expand"}
TEMPLATE_NAMES = {"asm", "scripts", "input", "output", "params", "expand"}
USER_AGENT = "pipeline/1.0"  # some CDNs reject urllib's default User-Agent


class RefError(Exception):
    """A catalog or spec problem; reported before any work starts."""


class Unavailable(Exception):
    """A catalog field the chosen assembly does not define."""


class Paths(list):
    """Paths that format as a space-separated string but can be indexed."""

    def __str__(self):
        return " ".join(self)

    def __format__(self, spec):
        return format(str(self), spec)


# -----------------------------------------------------------------------------
# References engine: loading
# -----------------------------------------------------------------------------
def _load_yaml(path):
    try:
        with open(path) as fh:
            return yaml.safe_load(fh) or {}
    except FileNotFoundError:
        raise RefError(f"file not found: {path}")
    except yaml.YAMLError as exc:
        raise RefError(f"{path}: invalid YAML: {exc}")


def load_catalog(path):
    """{assembly: {field: value}} with `defaults` merged under every entry."""
    data = _load_yaml(path)
    defaults = data.get("defaults") or {}
    assemblies = data.get("assemblies") or {}
    if not assemblies:
        raise RefError(f"{path}: no `assemblies:`")
    out = {}
    for name, fields in assemblies.items():
        if not isinstance(fields, dict):
            raise RefError(f"{path}: assembly {name!r} must be a mapping")
        out[str(name)] = {**defaults, **fields}
    return out


def load_port(name, pdir, images):
    """A port descriptor for plan_refs(): its spec, config and image paths.

    `images` maps each containers.yaml key to an absolute image path."""
    spec_path = os.path.join(pdir, "config", "references.yaml")
    spec = _load_yaml(spec_path)
    cfg_path = os.path.join(pdir, "config", "config.yaml")
    config = _load_yaml(cfg_path) if os.path.isfile(cfg_path) else {}
    return {"name": name, "dir": pdir, "spec_path": spec_path,
            "files": spec.get("files") or {}, "config": config,
            "images": images, "scripts": os.path.join(pdir, "scripts")}


def has_spec(pdir):
    return os.path.isfile(os.path.join(pdir, "config", "references.yaml"))


# -----------------------------------------------------------------------------
# References engine: planning
# -----------------------------------------------------------------------------
class Job:
    """One fetch or one build command; identified by its outputs."""

    def __init__(self, port, entry, kind, outputs, binding):
        self.port = port            # port descriptor (first port that needs it)
        self.entry = entry          # entry id in that port's spec
        self.kind = kind            # "fetch" | "build"
        self.outputs = outputs      # paths relative to the assembly folder
        self.binding = binding      # foreach variables
        self.source = None          # fetch: catalog source
        self.inputs = {}            # build: from-id -> relative paths
        self.image = None           # build: absolute image path
        self.container = None       # build: containers.yaml key
        self.shell = None           # build: shell template
        self.params = {}
        self.expand = {}
        self.unavailable = None     # reason, if it cannot be made
        self.pipelines = [port["name"]]

    @property
    def name(self):
        tag = ",".join(f"{k}={v}" for k, v in self.binding.items())
        return f"{self.port['name']}:{self.entry}" + (f"[{tag}]" if tag else "")


class Plan:
    def __init__(self, assembly, fields, repo):
        self.assembly = assembly
        self.fields = fields        # the catalog entry
        self.repo = repo
        self.jobs = []              # topological order
        self.producer = {}          # relative path -> Job
        self.file_fields = set()    # catalog fields used as files

    def pending(self, adir, force=False):
        """Jobs that execute_refs() would run (ignoring failures on the way)."""
        return [j for j in self.jobs if not j.unavailable
                and (force or not all(os.path.exists(os.path.join(adir, o)) for o in j.outputs))]

    def values(self):
        """Catalog fields that are plain values, not files."""
        return {k: v for k, v in self.fields.items()
                if k not in self.file_fields and k != "description" and not _is_source(v)}


def _dotted(data, path, what):
    cur = data
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            raise KeyError(f"{what} has no {path!r}")
        cur = cur[part]
    return cur


def _lookup(value, port, fields, binding):
    """Resolve a `config:` / `catalog:` reference; other values are literal."""
    if not isinstance(value, str):
        return value
    for prefix in ("config:", "catalog:"):
        if value.startswith(prefix):
            ref = value[len(prefix):].format(**binding)
            if prefix == "config:":
                try:
                    return _dotted(port["config"], ref, "config.yaml")
                except KeyError as exc:
                    raise RefError(f"{port['spec_path']}: {exc.args[0]}")
            if ref not in fields:
                raise Unavailable(ref)
            return fields[ref]
    return value


def _values(spec_value, port, fields, what):
    """A foreach/expand value list: a literal list, or a config:/catalog:
    reference to a list (its items) or a mapping (its keys)."""
    try:
        val = _lookup(spec_value, port, fields, {})
    except Unavailable as u:
        raise RefError(f"{port['spec_path']}: {what} reads catalog field {u.args[0]!r}, "
                       f"which this assembly lacks")
    if isinstance(val, dict):
        return list(val)
    if isinstance(val, (list, tuple)):
        return list(val)
    raise RefError(f"{port['spec_path']}: {what} must name a list or mapping, got {val!r}")


def _template_fields(tpl):
    """(top-level name, first index) of each replacement field in a template."""
    out = []
    for _, field, _, _ in string.Formatter().parse(tpl):
        if field is None:
            continue
        m = re.match(r"([A-Za-z_][A-Za-z0-9_]*)(?:\[([^\]]*)\])?", field)
        if not m:
            out.append((field, None))
        else:
            out.append((m.group(1), m.group(2)))
    return out


def _validate_entry(port, eid, e, ids):
    where = f"{port['spec_path']}: entry {eid!r}"
    if not isinstance(e, dict):
        raise RefError(f"{where} must be a mapping")
    if "out" not in e:
        raise RefError(f"{where}: missing `out`")
    if "catalog" in e:
        extra = set(e) - FETCH_KEYS
        if extra:
            raise RefError(f"{where}: a fetch entry (catalog:) takes only `catalog` and "
                           f"`out`; unexpected {', '.join(sorted(extra))}")
        if not isinstance(e["out"], str):
            raise RefError(f"{where}: a fetch entry has a single `out` path")
        return
    extra = set(e) - BUILD_KEYS
    if extra:
        raise RefError(f"{where}: unknown field(s) {', '.join(sorted(extra))}")
    for key in ("container", "shell"):
        if not e.get(key):
            raise RefError(f"{where}: a build entry needs `{key}` (or `catalog` for a fetch)")
    deps = e.get("from") or []
    if not isinstance(deps, list):
        raise RefError(f"{where}: `from` must be a list of entry ids")
    for d in deps:
        if d not in ids:
            raise RefError(f"{where}: `from` names unknown entry {d!r}")
    if e["container"] not in port["images"]:
        raise RefError(f"{where}: container {e['container']!r} is not a key of "
                       f"containers.yaml `containers:`")
    for key in ("foreach", "expand"):
        if key in e and not isinstance(e[key], dict):
            raise RefError(f"{where}: `{key}` must map variable names to value lists")
    if "params" in e and not isinstance(e["params"], dict):
        raise RefError(f"{where}: `params` must be a mapping")
    allowed = TEMPLATE_NAMES | set(e.get("foreach") or {})
    for name, index in _template_fields(e["shell"]):
        if name not in allowed:
            raise RefError(f"{where}: unknown placeholder {{{name}}} in shell "
                           f"(literal braces are written {{{{ }}}})")
        if name == "input" and index not in deps:
            raise RefError(f"{where}: {{input[{index}]}} is not listed in `from`")
        if name == "params" and index not in (e.get("params") or {}):
            raise RefError(f"{where}: {{params[{index}]}} is not defined in `params`")
        if name == "expand" and index not in (e.get("expand") or {}):
            raise RefError(f"{where}: {{expand[{index}]}} is not defined in `expand`")


def _entry_order(port):
    """Entry ids in dependency order; raises on a cycle."""
    files = port["files"]
    order, state = [], {}

    def visit(eid, stack):
        if state.get(eid) == "done":
            return
        if state.get(eid) == "active":
            cycle = " -> ".join(stack[stack.index(eid):] + [eid])
            raise RefError(f"{port['spec_path']}: dependency cycle {cycle}")
        state[eid] = "active"
        for d in files[eid].get("from") or []:
            visit(d, stack + [eid])
        state[eid] = "done"
        order.append(eid)

    for eid in files:
        visit(eid, [])
    return order


def _format(tpl, ns, where):
    try:
        return tpl.format(**ns)
    except (KeyError, IndexError, ValueError) as exc:
        raise RefError(f"{where}: cannot fill {tpl!r}: {exc}")


def plan_refs(catalog, assembly, ports, repo):
    """Expand the ports' entries into jobs for one assembly."""
    if assembly not in catalog:
        raise RefError(f"unknown assembly {assembly!r}; catalog has: {', '.join(catalog)}")
    fields = catalog[assembly]
    p = Plan(assembly, fields, repo)
    for port in ports:
        files = port["files"]
        if not isinstance(files, dict) or not files:
            raise RefError(f"{port['spec_path']}: no `files:`")
        for eid, e in files.items():
            _validate_entry(port, eid, e, files)
        entry_jobs = {}
        for eid in _entry_order(port):
            e = files[eid]
            where = f"{port['spec_path']}: entry {eid!r}"
            if "catalog" in e:
                p.file_fields.add(e["catalog"])
                bindings = [{}]
            else:
                fe = e.get("foreach") or {}
                lists = [_values(v, port, fields, f"foreach {k!r}") for k, v in fe.items()]
                bindings = [dict(zip(fe, combo)) for combo in itertools.product(*lists)]
            jobs = []
            for binding in bindings:
                jobs.append(_make_job(p, port, eid, e, binding, entry_jobs, where))
            entry_jobs[eid] = jobs
    return p


def _make_job(p, port, eid, e, binding, entry_jobs, where):
    fields, asm = p.fields, p.assembly
    unavailable = None
    expand = {}
    if "catalog" in e:
        outs = [_format(e["out"], {"asm": asm}, where)]
    else:
        ex = e.get("expand") or {}
        expand = {k: [str(x) for x in _values(v, port, fields, f"expand {k!r}")]
                  for k, v in ex.items()}
        out_tpls = e["out"] if isinstance(e["out"], list) else [e["out"]]
        outs = []
        combos = [dict(zip(expand, c)) for c in itertools.product(*expand.values())] or [{}]
        for tpl in out_tpls:
            names = {n for n, _ in _template_fields(tpl)}
            for combo in (combos if names & set(expand) else [{}]):
                path = _format(tpl, {"asm": asm, **binding, **combo}, where)
                if path not in outs:
                    outs.append(path)
    for o in outs:
        if os.path.isabs(o) or ".." in o.split("/"):
            raise RefError(f"{where}: output {o!r} must be relative to the assembly folder")

    # Same outputs as an earlier job (another port, usually): reuse it.
    known = [p.producer[o] for o in outs if o in p.producer]
    if known:
        first = known[0]
        if len(known) != len(outs) or any(k is not first for k in known) \
                or set(first.outputs) != set(outs):
            raise RefError(f"{where}: outputs {outs} overlap a different set made by "
                           f"{first.name} ({first.outputs})")
        if port["name"] not in first.pipelines:
            first.pipelines.append(port["name"])
        return first

    job = Job(port, eid, "fetch" if "catalog" in e else "build", outs, binding)
    if job.kind == "fetch":
        field = e["catalog"]
        if field in fields:
            job.source = str(fields[field])
        else:
            unavailable = f"catalog field {field!r} missing for {asm}"
    else:
        job.container = e["container"]
        job.image = port["images"][e["container"]]
        job.shell = e["shell"]
        job.expand = {k: Paths(v) for k, v in expand.items()}
        for k, v in (e.get("params") or {}).items():
            try:
                job.params[k] = _lookup(v, port, fields, binding)
            except Unavailable as u:
                unavailable = unavailable or f"catalog field {u.args[0]!r} missing for {asm}"
        for d in e.get("from") or []:
            dep_outs = [o for dj in entry_jobs[d] for o in dj.outputs]
            job.inputs[d] = dep_outs
            for dj in entry_jobs[d]:
                if dj.unavailable and not unavailable:
                    unavailable = f"needs {dj.name} ({dj.unavailable})"
    job.unavailable = unavailable
    p.jobs.append(job)
    for o in outs:
        p.producer[o] = job
    return job


# -----------------------------------------------------------------------------
# References engine: execution
# -----------------------------------------------------------------------------
def _is_source(value):
    """A file source (a field no selected port uses), not a plain value."""
    return isinstance(value, str) and (_is_url(value) or value.startswith(("repo:", "file://"))
                                       or os.path.isabs(value))


def _is_url(src):
    return re.match(r"^(https?|ftp)://", src) is not None


def _local_source(src, repo):
    if src.startswith("repo:"):
        return os.path.join(repo, src[len("repo:"):])
    if src.startswith("file://"):
        return src[len("file://"):]
    if os.path.isabs(src):
        return src
    raise RefError(f"source {src!r} is not a URL, repo:, file:// or absolute path")


def fetch_ref(src, dest, repo):
    """Download or copy `src` to `dest`, gunzipping a .gz source into a plain
    destination. Written to dest.part first, so a failure leaves no dest."""
    src_path = urllib.parse.urlparse(src).path if _is_url(src) else src
    gunzip = src_path.endswith(".gz") and not dest.endswith(".gz")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    part = dest + ".part"
    try:
        if _is_url(src):
            req = urllib.request.Request(src, headers={"User-Agent": USER_AGENT})
            stream = urllib.request.urlopen(req, timeout=60)
        else:
            stream = open(_local_source(src, repo), "rb")
        with stream as raw, open(part, "wb") as out:
            fh = gzip.GzipFile(fileobj=raw) if gunzip else raw
            shutil.copyfileobj(fh, out, length=1 << 20)
        os.replace(part, dest)
    except BaseException:
        if os.path.exists(part):
            os.remove(part)
        raise


def local_runner(image, cmd, cwd, binds):
    """Runs the command with the host's bash (tests; no container)."""
    return subprocess.call(["bash", "-euo", "pipefail", "-c", cmd], cwd=cwd)


def _safe(name):
    return re.sub(r"[^A-Za-z0-9_.=-]+", "_", name)


def build_ref(job, adir, runner):
    """Run a build job in <adir>/.tmp/<job>/; move its outputs in on success.
    Returns None, or an error message (the temp dir is then kept)."""
    tmp = os.path.join(adir, ".tmp", _safe(job.name))
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)
    tmp_outs = Paths(os.path.join(tmp, o) for o in job.outputs)
    for o in tmp_outs:
        os.makedirs(os.path.dirname(o), exist_ok=True)
    ns = {"asm": _assembly_of(adir), "scripts": job.port["scripts"],
          "input": {d: Paths(os.path.join(adir, o) for o in outs)
                    for d, outs in job.inputs.items()},
          "output": tmp_outs, "params": job.params, "expand": job.expand, **job.binding}
    cmd = job.shell.format(**ns)
    rc = runner(job.image, cmd, tmp, [adir, job.port["scripts"]])
    if rc != 0:
        return f"exit status {rc}; command and partial outputs kept in {tmp}"
    missing = [o for o, t in zip(job.outputs, tmp_outs) if not os.path.exists(t)]
    if missing:
        return f"command succeeded but did not write {', '.join(missing)} (kept {tmp})"
    for o, t in zip(job.outputs, tmp_outs):
        final = os.path.join(adir, o)
        os.makedirs(os.path.dirname(final), exist_ok=True)
        if os.path.isdir(final) and not os.path.islink(final):
            shutil.rmtree(final)
        os.replace(t, final)
    shutil.rmtree(tmp, ignore_errors=True)
    return None


def _assembly_of(adir):
    return os.path.basename(os.path.normpath(adir))


def execute_refs(p, refdir, runner, force=False, log=print):
    """Run the plan into <refdir>/<assembly>/. Returns {job name: status}."""
    adir = os.path.join(os.path.abspath(refdir), p.assembly)
    os.makedirs(adir, exist_ok=True)
    status = {}
    by_job = {}
    n = len(p.jobs)
    for i, job in enumerate(p.jobs, 1):
        head = f"[{i}/{n}] {job.name}"
        if job.unavailable:
            st = "unavailable"
            log(f"{head}: skipped, {job.unavailable}")
        elif any(by_job.get(id(p.producer[o])) in ("failed", "blocked")
                 for outs in job.inputs.values() for o in outs):
            st = "blocked"
            log(f"{head}: blocked by a failed input")
        elif not force and all(os.path.exists(os.path.join(adir, o)) for o in job.outputs):
            st = "existing"
            log(f"{head}: exists")
        elif job.kind == "fetch":
            dest = os.path.join(adir, job.outputs[0])
            log(f"{head}: {job.source} -> {job.outputs[0]}")
            sys.stdout.flush()
            try:
                fetch_ref(job.source, dest, p.repo)
                st = "downloaded"
            except Exception as exc:  # network, HTTP, disk, bad source
                st = "failed"
                log(f"    failed: {exc}")
        else:
            log(f"{head}: {', '.join(job.outputs)} ({job.container})")
            sys.stdout.flush()
            err = build_ref(job, adir, runner)
            st = "failed" if err else "built"
            if err:
                log(f"    failed: {err}")
        by_job[id(job)] = st
        status[job.name] = st
    write_readme(p, adir)
    try:
        os.rmdir(os.path.join(adir, ".tmp"))  # only when empty (no failed build kept)
    except OSError:
        pass
    return status


def write_readme(p, adir):
    """<adir>/README.tsv: every file present and every catalog value."""
    rows = [("kind", "name", "path", "pipelines", "source")]
    for job in p.jobs:
        for o in job.outputs:
            path = os.path.join(adir, o)
            if not os.path.exists(path):
                continue
            if job.kind == "fetch":
                src = job.source
            else:
                src = f"{job.port['name']} {job.container} ({os.path.basename(job.image)})"
            rows.append(("file", job.entry, path, ",".join(job.pipelines), src))
    for k, v in p.values().items():
        rows.append(("value", k, str(v), "", "catalog"))
    with open(os.path.join(adir, "README.tsv"), "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")


def refs_summary(status):
    counts = {}
    for st in status.values():
        counts[st] = counts.get(st, 0) + 1
    order = ["downloaded", "built", "existing", "unavailable", "failed", "blocked"]
    return ", ".join(f"{counts[s]} {s}" for s in order if s in counts)


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
        "references", help="download and build reference files",
        description="Download the base reference files of one assembly and build the "
                    "light files the ports derive from them (no aligner indexes) into "
                    "PATH/<assembly>/. Without --pipeline, the files of every port.")
    p.add_argument("-p", "--path", required=True, help="references directory, e.g. ~/references")
    p.add_argument("-a", "--assembly", help="catalog assembly (menu if omitted)")
    p.add_argument("--pipeline", action="append", default=[],
                   help="only this pipeline's files (repeatable)")
    p.add_argument("-n", "--dry-run", action="store_true",
                   help="only print what would be downloaded or built")
    p.add_argument("-f", "--force", action="store_true",
                   help="download/build again even if the files exist")
    p.set_defaults(func=cmd_references)

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
