# Every external tool pinned, GLAPD built from source at a pinned commit.
FROM mambaorg/micromamba:1.5.10

ARG GLAPD_COMMIT=0fcce4d39ccf6d4ff4e0b36045450a81ec0cf0dd

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates gcc g++ make \
 && rm -rf /var/lib/apt/lists/*

# GLAPD is GPL-2.0 and is invoked as a separate process; its source is never
# vendored into or patched by LAMPpy.
RUN git clone https://github.com/jiqingxiaoxi/GLAPD.git /opt/GLAPD \
 && cd /opt/GLAPD \
 && git checkout "$GLAPD_COMMIT" \
 && git rev-parse HEAD > /opt/GLAPD/COMMIT \
 && make \
 && rm -rf /opt/GLAPD/bowtie

COPY . /opt/lamppy
RUN pip install --no-cache-dir /opt/lamppy && rm -rf /opt/lamppy/.git

ENV GLAPD=/opt/GLAPD
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "lamppy"]
CMD ["--help"]
