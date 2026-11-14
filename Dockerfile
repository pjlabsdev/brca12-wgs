# brca12-wgs container: every tool of the pipeline with the versions pinned in envs/*.yml.
# Manta 1.6 needs Python 2.7 and lives in its own environment (/opt/conda/envs/manta).
# The Ensembl VEP cache, the reference genome and ClinVar are NOT included; mount them at run time.
#   docker build -t brca12-wgs:1.1.0 .
#   nextflow run main.nf -profile docker ...
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="brca12-wgs" \
      org.opencontainers.image.description="BRCA1/BRCA2 germline analysis from WGS, exome and panel data (GRCh38). Research use only." \
      org.opencontainers.image.source="https://github.com/pjlabsdev/brca12-wgs" \
      org.opencontainers.image.licenses="MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/brca12.yml envs/manta.yml /tmp/envs/
RUN micromamba install -y -n base -f /tmp/envs/brca12.yml \
 && micromamba install -y -n base -c conda-forge procps-ng \
 && micromamba create -y -n manta -f /tmp/envs/manta.yml \
 && micromamba clean -a -y \
 && rm -rf /tmp/envs

# The manta environment comes last on PATH: it only adds python2, which Manta's scripts call.
ENV PATH=/opt/conda/bin:$PATH:/opt/conda/envs/manta/bin \
    MANTA_CONFIG=/opt/conda/envs/manta/bin/configManta.py

# Smoke test at build time: every executable answers (outputs are captured first, so that a tool
# that exits non-zero on --help or a closed pipe does not hide the check).
RUN set -e; \
    v() { out=$("$@" 2>&1 || true); printf '%s\n' "$out" | sed -n 1p; [ -n "$out" ]; }; \
    out=$(bwa 2>&1 || true); printf '%s\n' "$out" | grep Version; \
    v samtools --version; v bcftools --version; v mosdepth --version; v cnvkit.py version; \
    out=$(gatk --version 2>&1 || true); printf '%s\n' "$out" | grep "Genome Analysis Toolkit"; \
    out=$(vep --help 2>&1 || true); printf '%s\n' "$out" | grep ensembl-vep; \
    v /opt/conda/envs/manta/bin/configManta.py --version; v ps --version
