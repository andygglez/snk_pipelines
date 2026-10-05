"""Software versions of a run, as nf-core reports them.

Used from the Snakefile's `software_versions` rule (run: block, host Python):

  collect()     for every nf-core process in config/versions.yaml whose port
                rules are part of this run's DAG (and whose `when` holds), run
                the module's version commands in the image of that rule
  write_topic() softwareVersionsToYAML + collectFile(sort: true):
                pipeline_info/nf_core_<name>_software_mqc_versions.yml
  write_dump()  CUSTOM_DUMPSOFTWAREVERSIONS (dumpsoftwareversions.py template):
                software_versions.yml + software_versions_mqc.yml

The Workflow section names Snakemake (and its version) where nf-core names
Nextflow.
"""

import platform
import re
import shutil
import subprocess
from collections import defaultdict
from textwrap import dedent

import yaml


def _run(img, cmd, cache):
    key = (img, cmd)
    if key not in cache:
        if img and shutil.which("apptainer"):
            argv = ["apptainer", "exec", img, "bash", "-c", cmd]
        elif img and shutil.which("singularity"):
            argv = ["singularity", "exec", img, "bash", "-c", cmd]
        else:
            argv = ["bash", "-c", cmd]
        try:
            res = subprocess.run(argv, capture_output=True, text=True, timeout=300)
            cache[key] = res.stdout.strip()
        except (OSError, subprocess.TimeoutExpired):
            cache[key] = ""
    return cache[key]


class _DagJob:
    def __init__(self, rule, wildcards, img):
        self.rule = type("Rule", (), {"name": rule})()
        self.wildcards_dict = wildcards
        self.container_img_path = img


def dag_jobs(workflow, files):
    """The jobs that made `files` and everything upstream of them.

    A `run:` job is executed by a separate Snakemake process whose
    workflow.dag holds only that job, so the jobs are traced through the
    metadata record Snakemake keeps for every output (producing rule and its
    inputs); wildcards come back from the rule's output patterns."""
    read = workflow.persistence._read_record
    seen, queue, jobs = set(), [str(f) for f in files], {}
    while queue:
        path = queue.pop()
        if path in seen:
            continue
        seen.add(path)
        record = read(path)
        if record is None or not record.rule:
            continue  # a source file
        rule = workflow.get_rule(record.rule)
        wildcards = {k: str(v) for k, v in dict(rule.get_wildcards(path)).items()}
        key = (record.rule, tuple(sorted(wildcards.items())))
        if key not in jobs:
            img = rule.container_img
            jobs[key] = _DagJob(record.rule, wildcards, img if isinstance(img, str) else None)
        queue.extend(record.input or [])
    return list(jobs.values())


def collect(table_path, jobs, namespace, yaml_numbers=True):
    """{process: {tool: version}} for the processes this DAG runs.

    yaml_numbers: render file-style versions as nf-core's SnakeYAML round
    trip does (softwareVersionsToYAML); False keeps the strings
    (CUSTOM_DUMPSOFTWAREVERSIONS reads versions.yml with yaml.BaseLoader)."""
    with open(table_path) as fh:
        table = yaml.safe_load(fh) or {}
    by_rule = defaultdict(list)
    for job in jobs:
        by_rule[job.rule.name].append(job)
    cache = {}
    versions = {}
    for process, entry in table.items():
        matched = [j for name, js in by_rule.items() for j in js
                   if any(re.fullmatch(p, name) for p in entry["rules"])]
        if not matched and entry.get("downstream"):
            # runs after this rule (MultiQC): present when its `when` holds
            rule = namespace["workflow"].get_rule(entry["rules"][0])
            img = rule.container_img
            matched = [_DagJob(rule.name, {}, img if isinstance(img, str) else None)]
        if not matched:
            continue
        if entry.get("when"):
            def any_wc(name, *values, _jobs=matched):
                return any(str(j.wildcards_dict.get(name, "")) in values for j in _jobs)

            env = dict(namespace)
            env.update(jobs=matched, any_wc=any_wc)
            if not eval(entry["when"], env):  # noqa: S307 (expressions from config/versions.yaml)
                continue
        try:
            img = matched[0].container_img_path
        except Exception:  # noqa: BLE001 (no container deployment)
            img = None
        found = {tool: _run(img, cmd, cache) for tool, cmd in entry["tools"]}
        if entry.get("style", "file") == "topic":
            # topic channels: one (process, tool, version) tuple per tool, grouped sorted
            versions[process] = dict(sorted(found.items()))
        else:
            # versions.yml: SnakeYAML load + dumpAsMap keeps the module's order
            # and turns numeric-looking versions into numbers (1.20 -> 1.2)
            versions[process] = {tool: _yaml_scalar(v) if yaml_numbers else v for tool, v in found.items()}
    return versions


def _yaml_scalar(value):
    """A version string as SnakeYAML load + dump would print it."""
    try:
        parsed = yaml.safe_load(value)
    except yaml.YAMLError:
        return value
    if isinstance(parsed, bool) or parsed is None:
        return value
    if isinstance(parsed, (int, float)):
        return repr(parsed)
    return value


def _snakemake_version():
    try:
        import snakemake

        return snakemake.__version__
    except Exception:  # noqa: BLE001
        return "unknown"


def version_chunks(versions):
    """One 'PROCESS:\n  tool: version' chunk per process."""
    return [f"{process}:\n" + "\n".join(f"  {t}: {v}" for t, v in tools.items())
            for process, tools in versions.items()]


def workflow_chunk(name, version):
    return f"Workflow:\n    {name}: {version}\n    Snakemake: {_snakemake_version()}"


def collect_file(chunks):
    """collectFile(sort: true, newLine: true) of the chunks."""
    return "".join(c + "\n" for c in sorted(set(chunks)))


def styles(table_path):
    with open(table_path) as fh:
        return {p: e.get("style", "file") for p, e in (yaml.safe_load(fh) or {}).items()}


def write_topic(versions, name, version, path):
    """nf_core_<name>_software_mqc_versions.yml (chunks sorted, one newline each)."""
    with open(path, "w") as fh:
        fh.write(collect_file(version_chunks(versions) + [workflow_chunk(name, version)]))


def _make_versions_html(versions):
    html = [dedent(
        """\
        <style>
        #nf-core-versions tbody:nth-child(even) {
            background-color: #f2f2f2;
        }
        </style>
        <table class="table" style="width:100%" id="nf-core-versions">
            <thead>
                <tr>
                    <th> Process Name </th>
                    <th> Software </th>
                    <th> Version  </th>
                </tr>
            </thead>
        """
    )]
    for process, tmp_versions in sorted(versions.items()):
        html.append("<tbody>")
        for i, (tool, ver) in enumerate(sorted(tmp_versions.items())):
            html.append(dedent(
                f"""\
                <tr>
                    <td><samp>{process if (i == 0) else ''}</samp></td>
                    <td><samp>{tool}</samp></td>
                    <td><samp>{ver}</samp></td>
                </tr>
                """
            ))
        html.append("</tbody>")
    html.append("</table>")
    return "\n".join(html)


def write_dump(versions, name, version, yml_path, mqc_path, img=None):
    """CUSTOM_DUMPSOFTWAREVERSIONS outputs (versions as strings, BaseLoader);
    its own python/yaml versions come from its image (img)."""
    by_module = dict(versions)
    cache = {}
    py = _run(img, "python -c 'import platform; print(platform.python_version())'", cache) if img else ""
    ym = _run(img, "python -c 'import yaml; print(yaml.__version__)'", cache) if img else ""
    by_module["CUSTOM_DUMPSOFTWAREVERSIONS"] = {"python": py or platform.python_version(), "yaml": ym or yaml.__version__}
    by_module["Workflow"] = {"Snakemake": _snakemake_version(), name: version}
    mqc = {
        "id": "software_versions",
        "section_name": f"{name} Software Versions",
        "section_href": f"https://github.com/{name}",
        "plot_type": "html",
        "description": "are collected at run time from the software output.",
        "data": _make_versions_html(by_module),
    }
    with open(yml_path, "w") as fh:
        yaml.dump(by_module, fh, default_flow_style=False)
    with open(mqc_path, "w") as fh:
        yaml.dump(mqc, fh, default_flow_style=False)
