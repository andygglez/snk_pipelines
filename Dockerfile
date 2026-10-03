# Runs the `pipeline` command (pipeline.py) with Snakemake installed by uv and
# Apptainer for the per-rule Singularity images.
#
#   docker build -t pipelines .
#   docker run --rm -it --privileged --user "$(id -u):$(id -g)" \
#       -v "$HOME/pipeline-images":/containers \
#       -v "$PWD":/work \
#       pipelines run rnaseq -c rnaseq.config.yaml -t 8
#
# --user keeps the outputs owned by you instead of root.
# Images are not baked in: each snk-*/containers/ is a symlink to
# /containers/snk-*/, so mount one volume (or host directory) at /containers and
# fill it once with `pipeline pull <name>`. --privileged lets Apptainer run
# inside Docker (setuid starter + mounts of the images).
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim

ARG APPTAINER_VERSION=1.3.5

# Apptainer from the upstream setuid .deb (same version as Rorqual's module):
# setuid mode works for non-root users even where the host blocks unprivileged
# user namespaces (e.g. Ubuntu's AppArmor restriction). procps
# provides `ps`, which Snakemake uses.
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl procps \
 && base=https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VERSION} \
 && curl -fsSL -o /tmp/apptainer.deb "$base/apptainer_${APPTAINER_VERSION}_amd64.deb" \
 && curl -fsSL -o /tmp/apptainer-suid.deb "$base/apptainer-suid_${APPTAINER_VERSION}_amd64.deb" \
 && apt-get install -y --no-install-recommends /tmp/apptainer.deb /tmp/apptainer-suid.deb \
 && rm -rf /tmp/apptainer*.deb /var/lib/apt/lists/*

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    PATH=/opt/venv/bin:$PATH

WORKDIR /opt/pipelines

# Dependencies first, so editing the ports doesn't reinstall Snakemake.
COPY pyproject.toml uv.lock README.md ./
RUN uv sync --frozen --no-install-project --no-dev

# The ports, then the project itself (editable: pipeline.py must stay next to snk-*/).
COPY . .
RUN uv sync --frozen --no-dev \
 && for d in snk-*/; do \
      d=${d%/}; rm -rf "$d/containers"; mkdir -p "/containers/$d"; \
      ln -s "/containers/$d" "$d/containers"; \
    done

# --user $(id -u):$(id -g) runs as an arbitrary uid: Apptainer needs a passwd
# entry for it, so let the entrypoint add one (container-local files only).
RUN chmod 666 /etc/passwd /etc/group

COPY --chmod=755 <<'EOF' /usr/local/bin/entrypoint.sh
#!/bin/sh
# Non-root (docker run --user): register the uid/gid and use a writable HOME.
if [ "$(id -u)" != 0 ]; then
    export HOME=/tmp/home
    mkdir -p "$HOME"
    getent group "$(id -g)" >/dev/null || echo "user:x:$(id -g):" >> /etc/group
    getent passwd "$(id -u)" >/dev/null \
        || echo "user:x:$(id -u):$(id -g)::$HOME:/bin/sh" >> /etc/passwd
fi
# A host directory mounted at /containers starts empty: recreate the per-port
# folders the snk-*/containers symlinks point to.
for d in /opt/pipelines/snk-*/; do mkdir -p "/containers/$(basename "$d")"; done
exec pipeline "$@"
EOF

WORKDIR /work
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["--help"]
