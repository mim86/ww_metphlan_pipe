# Wastewater MetaPhlAn + HoCoRT Snakemake pipeline (Docker)

This repository contains a Dockerised Snakemake pipeline for **Illumina paired-end metagenomic FASTQs** that performs:

Raw QC → adapter trimming → host (human) read removal → quality trimming → post-clean QC → MetaPhlAn profiling → merged abundance table.

It is designed so that everything can be run from a single Docker container, while keeping all inputs/outputs on your local disk via bind mounts.

## What’s included

The container environment includes Snakemake (v7.32.4), HoCoRT (v1.2.2), MetaPhlAn (v4.2.2), fastp (v1.0.1), FastQC (v0.12.1) and MultiQC (v1.33) (plus common mapping/utility tools).

Many paramters are controllable through a config.yaml file (threads, host index prefix, fastp options, MetaPhlAn DB/index/options).

---

## Folder layout

Expected (example):

```text
.
├── config.yaml
├── snakefile
├── raw_data/
│   └── RUN1 (current run)/
│       ├── sample1_R1.fastq.gz
│       ├── sample1_R2.fastq.gz
│       ├── sample2_R1.fastq.gz
│       └── sample2_R2.fastq.gz
├── human_genome/
│   ├── hg38.1.bt2
│   ├── hg38.2.bt2
│   ├── hg38.3.bt2
│   ├── hg38.4.bt2
│   ├── hg38.rev.1.bt2
│   └── hg38.rev.2.bt2
├── metaphlan_db (currently newest db available)
│   ├── mpa_latest
│   ├── mpa_vJan25_CHOCOPhlAnSGB_202503.1.bt2l
│   ├── mpa_vJan25_CHOCOPhlAnSGB_202503.2.bt2l
│   ├── mpa_vJan25_CHOCOPhlAnSGB_202503.3.bt2l
│   ├── mpa_vJan25_CHOCOPhlAnSGB_202503.4.bt2l
│   ├── mpa_vJan25_CHOCOPhlAnSGB_202503.nwk
│   ├── mpa_vJan25_CHOCOPhlAnSGB_202503.pkl
│   ├── mpa_vJan25_CHOCOPhlAnSGB_202503.rev.1.bt2l
│   ├── mpa_vJan25_CHOCOPhlAnSGB_202503.rev.2.bt2l
│   ├── mpa_vJan25_CHOCOPhlAnSGB_202503_VINFO.csv
│   └── mpa_vJan25_CHOCOPhlAnSGB_202503_VSG.fna
```

Your results will be written to:

```text
results/<RUN_NAME>/
```

---

## Quickstart (using the prebuilt Docker image)

### 1) Pull the Docker image

```bash
docker pull mimsto86/wastewater-metaphlan-hocort:v1
```

### 2) Install the MetaPhlAn database (one-time)

Create the DB directory:

```bash
mkdir -p metaphlan_db
```

Download/install the DB into `./metaphlan_db`:

```bash
docker run --rm -it \
  --user "$(id -u)":"$(id -g)" \
  -v "$(pwd)":/work -w /work \
  mimsto86/wastewater-metaphlan-hocort:v1 \
  metaphlan --install --db_dir /work/metaphlan_db
```

Notes:
- The `--user "$(id -u)":"$(id -g)"` is important so the database files are owned by you (and not root).
- The config expects the DB in `/work/metaphlan_db` inside the container, which corresponds to `./metaphlan_db` on your host.

### 3) Provide a Bowtie2 host (human) index

Config points to:

```yaml
host_index_prefix: "human_genome/hg38"
```

So the files must exist like:

```text
human_genome/hg38.1.bt2  ... hg38.rev.2.bt2
```

(These are Bowtie2 index files; the pipeline uses HoCoRT with Bowtie2 for filtering, although HoCoRT offers flexibility in terms which mapper could be used for filtering, please refer to [HoCoRT: host contamination removal tool](https://link.springer.com/article/10.1186/s12859-023-05492-w))

### 4) Put FASTQs in the expected naming format

The pipeline auto-detects samples by scanning:

```text
raw_data/<RUN>/*_R1.fastq.gz
```

For each `*_R1.fastq.gz`, it expects a matching `*_R2.fastq.gz` with the same sample prefix.

Example:
- `raw_data/RUN1/sample1_R1.fastq.gz`
- `raw_data/RUN1/sample1_R2.fastq.gz`

### 5) Run the pipeline

From the repo root (where `snakefile` and `config.yaml` are located):

```bash
docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$PWD":/work \
  -w /work \
  -e HOME=/work \
  -e XDG_CACHE_HOME=/work/.cache \
  mimsto86/wastewater-metaphlan-hocort:v1 \
  snakemake \
    --snakefile snakefile \
    --configfile config.yaml \
    --cores 8
```

---

## Configuration

All main settings are in `config.yaml`.

### Run selection

```yaml
current_run: "RUN1"
```

This tells the pipeline to look in `raw_data/RUN1/` and write results to `results/RUN1/`.

### Threads per tool

```yaml
threads:
  fastqc: 4
  fastp: 8
  hocort: 8
  metaphlan: 8
  multiqc: 2
```

### Host index prefix

```yaml
host_index_prefix: "human_genome/hg38"
```

### fastp behaviour

There are two fastp steps:

1) **Adapter-only trimming** (no quality/length filtering)  
Controlled by `fastp_adapter:`.

2) **Quality trimming + length filtering** on non-host reads  
Controlled by `fastp_clean:`.

Any `extra_args:` string is appended directly to fastp (so you can add custom flags).

### MetaPhlAn behaviour

```yaml
metaphlan:
  db_dir: "/work/metaphlan_db"
  index: "mpa_vJan25_CHOCOPhlAnSGB_202503"
  tax_lev: "a"
  analysis_type: "rel_ab"
  # optional cutoffs can be empty to use MetaPhlAn defaults:
  read_min_len:
  min_mapq_val:
  stat_q:
  perc_nonzero:
  # optional flags:
  skip_unclassified_estimation:
  subsampling_paired:
  subsampling_seed:
  extra_args:
```

If an optional value is blank, the pipeline does **not** add that flag, and MetaPhlAn will use its internal default.

---

## Pipeline steps and outputs

All outputs go to:

```text
results/<RUN>/
```

### 1) Raw QC (FastQC + MultiQC)

Directory:

```text
results/<RUN>/qc_raw/
```

Files:
- `*_R1_fastqc.html`, `*_R2_fastqc.html` (+ zipped FastQC outputs)
- `multiqc_raw.html` (+ `multiqc_raw_data/`)

Purpose:
- Baseline QC of the raw input FASTQs.

### 2) Adapter-only trimming (fastp)

Directory:

```text
results/<RUN>/trim_adapter/
```

Files per sample:
- `*_R1.trimA.fastq.gz`, `*_R2.trimA.fastq.gz`
- `*.trimA.fastp.html`, `*.trimA.fastp.json`

Purpose:
- Remove adapters (and only adapters) before host filtering.

### 3) Host removal (HoCoRT + Bowtie2)

Directory:

```text
results/<RUN>/host_removed/
```

Files per sample:
- `*_R1.nohost.fastq.gz`, `*_R2.nohost.fastq.gz`

Purpose:
- Map reads to the host index and keep **unmapped** reads (non-host).

### 4) Quality trimming / cleaning (fastp)

Directory:

```text
results/<RUN>/trim_clean/
```

Files per sample:
- `*_R1.clean.fastq.gz`, `*_R2.clean.fastq.gz`
- `*.clean.fastp.html`, `*.clean.fastp.json`

Purpose:
- Clean the non-host reads (quality trimming + length filtering).

Safety behaviour:
- If a sample ends up with **0 reads after host removal**, the pipeline skips fastp and propagates empty FASTQs while still producing minimal JSON/HTML reports (so downstream rules don’t crash).

### 5) Clean QC (FastQC + MultiQC)

Directory:

```text
results/<RUN>/qc_clean/
```

Files:
- `*_R1.clean_fastqc.html`, `*_R2.clean_fastqc.html`
- `multiqc_clean.html` (+ `multiqc_clean_data/`)

Purpose:
- QC after cleaning, so you can see what trimming/filtering did.

### 6) Taxonomic profiling (MetaPhlAn)

Directory:

```text
results/<RUN>/metaphlan/
```

Files per sample:
- `*_profile.txt` (MetaPhlAn profile output)
- `*.mapout.bz2` (MetaPhlAn bowtie2 mapping output)

Safety behaviour:
- If a sample has **0 cleaned reads**, the pipeline creates a stub `*_profile.txt` and an empty `*.mapout.bz2`.

### 7) Merge per-sample MetaPhlAn tables

File:

```text
results/<RUN>/metaphlan/merged_profiles.tsv
```

Purpose:
- A single merged abundance table across all samples.

---

## Building the Docker image yourself (optional)

If you want to rebuild the container locally (e.g., to modify versions or add tools), use the Docker recipe in `docker_file/` and the conda environment spec.

Example build:

```bash
cd docker_file
docker build -t my-wastewater-metaphlan-hocort:dev .
```

Then run exactly like the prebuilt image, replacing the image name:

```bash
docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$PWD":/work -w /work \
  -e HOME=/work -e XDG_CACHE_HOME=/work/.cache \
  my-wastewater-metaphlan-hocort:dev \
  snakemake --snakefile snakefile --configfile config.yaml --cores 8
```

---

## Troubleshooting

### MetaPhlAn DB not found / wrong index name
- Confirm `metaphlan_db/` exists on the host and contains the installed database files.
- Confirm your config matches the index you downloaded (example index in config is `mpa_vJan25_CHOCOPhlAnSGB_202503`).

### Permission denied when downloading DB
Use the provided `--user "$(id -u)":"$(id -g)"` and bind-mount the repo to `/work` so files are written with your UID/GID.

### “No samples found”
Ensure the raw FASTQs follow:

```text
raw_data/<RUN>/<sample>_R1.fastq.gz
raw_data/<RUN>/<sample>_R2.fastq.gz
```

…and that `current_run: "<RUN>"` points to the right folder.

### Snakemake is reusing old outputs unexpectedly
Try:

```bash
snakemake --snakefile snakefile --configfile config.yaml --cores 8 --rerun-incomplete
```

If you previously had a crash/interrupt and Snakemake locked the directory, run `--unlock`.
