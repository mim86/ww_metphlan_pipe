import os
import glob
from snakemake.io import glob_wildcards

# -----------------------------
# Load configuration
# -----------------------------
configfile: "config.yaml"

RUN = config["current_run"]
results_dir = f"results/{RUN}"

# Threads
THREADS = config.get("threads", {})
T_FASTQC    = THREADS.get("fastqc", 4)
T_FASTP     = THREADS.get("fastp", 8)
T_HOCORT    = THREADS.get("hocort", 8)
T_METAPHLAN = THREADS.get("metaphlan", 8)
T_MULTIQC   = THREADS.get("multiqc", 2)

# HoCoRT / Bowtie2 index prefix
HOST_INDEX_PREFIX = config["host_index_prefix"]

# Tool-specific configs
fastp_adapter_cfg = config.get("fastp_adapter", {})
fastp_clean_cfg   = config.get("fastp_clean", {})
metaphlan_cfg     = config.get("metaphlan", {})

# -----------------------------
# NEW: Accept multiple FASTQ extensions in raw_data
# -----------------------------
RAW_EXTS = ["fastq.gz", "fastq", "fq.gz", "fq"]

def get_raw_fastq(sample, read):
    """
    Return the existing raw FASTQ path for a given sample/read across:
    .fastq, .fastq.gz, .fq, .fq.gz

    Expected naming stays the same:
      {sample}_R1.<ext>
      {sample}_R2.<ext>
    """
    for ext in RAW_EXTS:
        p = f"raw_data/{RUN}/{sample}_R{read}.{ext}"
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"No raw FASTQ found for sample={sample} read=R{read} with extensions {RAW_EXTS} "
        f"in raw_data/{RUN}/"
    )

def _truthy(val):
    """
    Treat yes/true/True as True, and None/""/False/no as False.
    """
    if isinstance(val, bool):
        return val
    if val is None:
        return False
    s = str(val).strip().lower()
    return s in ("y", "yes", "true", "1")

def build_fastp_adapter_args():
    """
    Build optional fastp args for adapter-only trimming (trim_adapters).
    """
    parts = []

    if _truthy(fastp_adapter_cfg.get("detect_adapter_for_pe")):
        parts.append("--detect_adapter_for_pe")
    if _truthy(fastp_adapter_cfg.get("disable_quality_filtering")):
        parts.append("--disable_quality_filtering")
    if _truthy(fastp_adapter_cfg.get("disable_length_filtering")):
        parts.append("--disable_length_filtering")

    extra = fastp_adapter_cfg.get("extra_args")
    if extra not in (None, ""):
        parts.append(str(extra))

    return " ".join(parts)

def build_fastp_clean_args():
    """
    Build optional fastp args for quality trimming (trim_clean).
    """
    parts = []

    # boolean flags
    if _truthy(fastp_clean_cfg.get("detect_adapter_for_pe")):
        parts.append("--detect_adapter_for_pe")
    if _truthy(fastp_clean_cfg.get("cut_front")):
        parts.append("--cut_front")
    if _truthy(fastp_clean_cfg.get("cut_tail")):
        parts.append("--cut_tail")

    # numeric options (only add if user set them)
    val = fastp_clean_cfg.get("cut_tail_mean_quality")
    if val not in (None, ""):
        parts.append(f"--cut_tail_mean_quality {val}")

    val = fastp_clean_cfg.get("length_required")
    if val not in (None, ""):
        parts.append(f"--length_required {val}")

    extra = fastp_clean_cfg.get("extra_args")
    if extra not in (None, ""):
        parts.append(str(extra))

    return " ".join(parts)

def build_metaphlan_optional_args():
    """
    Build optional MetaPhlAn CLI args based on config.
    If a key is missing or empty, the flag is not added and
    MetaPhlAn's internal default is used.
    """
    cfg = metaphlan_cfg
    parts = []

    def add_if_present(key, flag):
        val = cfg.get(key)
        if val not in (None, "", False):
            parts.append(f"{flag} {val}")

    # numeric/string options
    add_if_present("read_min_len", "--read_min_len")
    add_if_present("min_mapq_val", "--min_mapq_val")
    add_if_present("stat_q", "--stat_q")
    add_if_present("perc_nonzero", "--perc_nonzero")

    # boolean flag
    if _truthy(cfg.get("skip_unclassified_estimation")):
        parts.append("--skip_unclassified_estimation")

    # subsampling options
    add_if_present("subsampling_paired", "--subsampling_paired")
    add_if_present("subsampling_seed", "--subsampling_seed")

    # free-text extras
    extra = cfg.get("extra_args")
    if extra not in (None, ""):
        parts.append(str(extra))

    return " ".join(parts)

METAPHLAN_DB_DIR = metaphlan_cfg["db_dir"]
METAPHLAN_INDEX  = metaphlan_cfg["index"]

# -----------------------------
# Helper: count reads in (possibly gzipped) FASTQ
# -----------------------------
def count_fastq_reads(path):
    """
    Count reads in a FASTQ (gzipped or plain).
    Returns number of reads (lines / 4).
    """
    import gzip
    if path.endswith(".gz"):
        opener = gzip.open
        mode = "rt"
    else:
        opener = open
        mode = "r"

    n_lines = 0
    with opener(path, mode, errors="ignore") as fh:
        for _ in fh:
            n_lines += 1
    return n_lines // 4

# -----------------------------
# Detect samples from raw R1 files (NOW supports fastq/fastq.gz/fq/fq.gz)
# -----------------------------
r1_files = []
for ext in RAW_EXTS:
    r1_files.extend(glob.glob(f"raw_data/{RUN}/*_R1.{ext}"))

SAMPLES = sorted({os.path.basename(p).split("_R1.")[0] for p in r1_files})

# -----------------------------
# Final targets
# -----------------------------
rule all:
    input:
        # MultiQC reports
        f"{results_dir}/qc_raw/multiqc_raw.html",
        f"{results_dir}/qc_clean/multiqc_clean.html",
        # Merged MetaPhlAn profiles
        f"{results_dir}/metaphlan/merged_profiles.tsv"


# =============================
# RAW QC
# =============================

rule fastqc_raw:
    """
    FastQC on raw reads (before any trimming).
    """
    input:
        r1 = lambda wc: get_raw_fastq(wc.sample, 1),
        r2 = lambda wc: get_raw_fastq(wc.sample, 2)
    output:
        html_r1 = f"{results_dir}/qc_raw/{{sample}}_R1_fastqc.html",
        html_r2 = f"{results_dir}/qc_raw/{{sample}}_R2_fastqc.html"
    threads: T_FASTQC
    shell:
        r"""
        mkdir -p {results_dir}/qc_raw
        fastqc -t {threads} -o {results_dir}/qc_raw {input.r1} {input.r2}
        """


rule multiqc_raw:
    """
    MultiQC summarising raw FastQC.
    """
    input:
        expand(f"{results_dir}/qc_raw/{{sample}}_R1_fastqc.html", sample=SAMPLES),
        expand(f"{results_dir}/qc_raw/{{sample}}_R2_fastqc.html", sample=SAMPLES)
    output:
        f"{results_dir}/qc_raw/multiqc_raw.html"
    threads: T_MULTIQC
    shell:
        r"""
        mkdir -p {results_dir}/qc_raw
        cd {results_dir}/qc_raw
        multiqc . -n multiqc_raw
        """


# =============================
# ADAPTER-ONLY TRIMMING (fastp)
# =============================

rule trim_adapters:
    """
    Adapter-only trimming with fastp (no quality / length filtering).
    """
    input:
        r1 = lambda wc: get_raw_fastq(wc.sample, 1),
        r2 = lambda wc: get_raw_fastq(wc.sample, 2)
    output:
        r1   = f"{results_dir}/trim_adapter/{{sample}}_R1.trimA.fastq.gz",
        r2   = f"{results_dir}/trim_adapter/{{sample}}_R2.trimA.fastq.gz",
        html = f"{results_dir}/trim_adapter/{{sample}}.trimA.fastp.html",
        json = f"{results_dir}/trim_adapter/{{sample}}.trimA.fastp.json"
    threads: T_FASTP
    params:
        fastp_args = build_fastp_adapter_args()
    shell:
        r"""
        mkdir -p {results_dir}/trim_adapter
        fastp \
          -i {input.r1} \
          -I {input.r2} \
          -o {output.r1} \
          -O {output.r2} \
          {params.fastp_args} \
          --thread {threads} \
          --html {output.html} \
          --json {output.json}
        """


# =============================
# HOST REMOVAL WITH HoCoRT
# =============================

rule hocort_host_removal:
    """
    Remove host reads using HoCoRT + Bowtie2, keeping unmapped (non-host) reads.
    """
    input:
        r1 = f"{results_dir}/trim_adapter/{{sample}}_R1.trimA.fastq.gz",
        r2 = f"{results_dir}/trim_adapter/{{sample}}_R2.trimA.fastq.gz"
    output:
        r1 = f"{results_dir}/host_removed/{{sample}}_R1.nohost.fastq.gz",
        r2 = f"{results_dir}/host_removed/{{sample}}_R2.nohost.fastq.gz"
    params:
        index = HOST_INDEX_PREFIX
    threads: T_HOCORT
    shell:
        r"""
        mkdir -p {results_dir}/host_removed
        hocort map bowtie2 \
          -x {params.index} \
          -i {input.r1} {input.r2} \
          -o {output.r1} {output.r2} \
          --filter true
        """


# =============================
# REAL QC: QUALITY TRIMMING
# =============================

rule trim_clean:
    """
    Quality trimming + length filtering on non-host reads (fastp).
    With safety: if both nohost FASTQs have 0 reads, skip fastp and just
    propagate empty files + simple reports.
    """
    input:
        r1 = f"{results_dir}/host_removed/{{sample}}_R1.nohost.fastq.gz",
        r2 = f"{results_dir}/host_removed/{{sample}}_R2.nohost.fastq.gz"
    output:
        r1   = f"{results_dir}/trim_clean/{{sample}}_R1.clean.fastq.gz",
        r2   = f"{results_dir}/trim_clean/{{sample}}_R2.clean.fastq.gz",
        html = f"{results_dir}/trim_clean/{{sample}}.clean.fastp.html",
        json = f"{results_dir}/trim_clean/{{sample}}.clean.fastp.json"
    threads: T_FASTP
    params:
        fastp_args = build_fastp_clean_args()
    run:
        from snakemake.shell import shell as sh

        n1 = count_fastq_reads(input.r1)
        n2 = count_fastq_reads(input.r2)

        sh(f"mkdir -p {results_dir}/trim_clean")

        if n1 == 0 and n2 == 0:
            # No reads after host removal – just propagate empties & write dummy reports
            sh(f"cp {input.r1} {output.r1}")
            sh(f"cp {input.r2} {output.r2}")
            # Minimal JSON/HTML so downstream steps don't choke
            with open(output.json, "w") as jf:
                jf.write('{"note": "no reads after host removal; fastp not run"}\n')
            with open(output.html, "w") as hf:
                hf.write(
                    f"<html><body><h1>No reads after host removal for sample {wildcards.sample}</h1></body></html>\n"
                )
        else:
            sh(
                f"fastp "
                f"-i {input.r1} "
                f"-I {input.r2} "
                f"-o {output.r1} "
                f"-O {output.r2} "
                f"{params.fastp_args} "
                f"--thread {threads} "
                f"--html {output.html} "
                f"--json {output.json}"
            )


# =============================
# QC AFTER CLEANING
# =============================

rule fastqc_clean:
    """
    FastQC on cleaned, non-host reads.
    Note: FastQC names outputs based on the full basename, so we expect *_R1.clean_fastqc.html.
    """
    input:
        r1 = f"{results_dir}/trim_clean/{{sample}}_R1.clean.fastq.gz",
        r2 = f"{results_dir}/trim_clean/{{sample}}_R2.clean.fastq.gz"
    output:
        html_r1 = f"{results_dir}/qc_clean/{{sample}}_R1.clean_fastqc.html",
        html_r2 = f"{results_dir}/qc_clean/{{sample}}_R2.clean_fastqc.html"
    threads: T_FASTQC
    shell:
        r"""
        mkdir -p {results_dir}/qc_clean
        fastqc -t {threads} -o {results_dir}/qc_clean {input.r1} {input.r2}
        """


rule multiqc_clean:
    """
    MultiQC summarising FastQC on cleaned reads.
    """
    input:
        expand(f"{results_dir}/qc_clean/{{sample}}_R1.clean_fastqc.html", sample=SAMPLES),
        expand(f"{results_dir}/qc_clean/{{sample}}_R2.clean_fastqc.html", sample=SAMPLES)
    output:
        f"{results_dir}/qc_clean/multiqc_clean.html"
    threads: T_MULTIQC
    shell:
        r"""
        mkdir -p {results_dir}/qc_clean
        cd {results_dir}/qc_clean
        multiqc . -n multiqc_clean
        """


# =============================
# METAPHLAN PROFILING PER SAMPLE
# =============================

rule metaphlan_profile:
    """
    Run MetaPhlAn on cleaned, non-host reads.
    With safety: if both cleaned FASTQs have 0 reads, skip MetaPhlAn and
    write a stub profile + empty mapout.
    """
    input:
        r1 = f"{results_dir}/trim_clean/{{sample}}_R1.clean.fastq.gz",
        r2 = f"{results_dir}/trim_clean/{{sample}}_R2.clean.fastq.gz"
    output:
        profile = f"{results_dir}/metaphlan/{{sample}}_profile.txt",
        bt2     = f"{results_dir}/metaphlan/{{sample}}.mapout.bz2"
    params:
        db_dir        = METAPHLAN_DB_DIR,
        index         = METAPHLAN_INDEX,
        tax_lev       = metaphlan_cfg.get("tax_lev", "a"),
        analysis_type = metaphlan_cfg.get("analysis_type", "rel_ab"),
        opt_args      = build_metaphlan_optional_args()
    threads: T_METAPHLAN
    run:
        from snakemake.shell import shell as sh

        n1 = count_fastq_reads(input.r1)
        n2 = count_fastq_reads(input.r2)

        sh(f"mkdir -p {results_dir}/metaphlan")

        if n1 == 0 and n2 == 0:
            # No reads after cleaning – don't run MetaPhlAn, create stub outputs
            sh(f"touch {output.bt2}")
            with open(output.profile, "w") as pf:
                pf.write("# MetaPhlAn profile stub: no reads after filtering; MetaPhlAn not run\n")
                pf.write("#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n")
        else:
            opt_args = params.opt_args or ""
            sh(
                f"metaphlan {input.r1},{input.r2} "
                f"--input_type fastq "
                f"--db_dir {params.db_dir} "
                f"-x {params.index} "
                f"--tax_lev {params.tax_lev} "
                f"-t {params.analysis_type} "
                f"{opt_args} "
                f"--mapout {output.bt2} "
                f"--nproc {threads} "
                f"-o {output.profile}"
            )


# =============================
# MERGE METAPHLAN TABLES
# =============================

rule metaphlan_merge:
    """
    Merge all per-sample MetaPhlAn profiles into a single abundance table.
    """
    input:
        expand(f"{results_dir}/metaphlan/{{sample}}_profile.txt", sample=SAMPLES)
    output:
        f"{results_dir}/metaphlan/merged_profiles.tsv"
    shell:
        r"""
        mkdir -p {results_dir}/metaphlan
        merge_metaphlan_tables.py {input} > {output}
        """

