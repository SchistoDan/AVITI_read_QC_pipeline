---
title: "SOP: Running the AVITI Read QC Pipeline"
author: "Dan Parsons @NHMUK"
date: "`14.04.2026`"
---

## Overview

This simple SOP describes how to configure and run the `aviti_read_qc_pipeline` (a Snakemake-based
pipeline for QC of raw, basecalled, and demultiplexed AVITI24 sequencing data). The pipeline
merges lane replicates, runs fastp trimming, falco (FastQC) reporting, and seqkit stats, and
aggregates everything into a fastp summary spreadsheet and MultiQC report. It can optionally
also align reads to a reference genome and report mapping statistics.

**Prerequisites:** conda must be installed and the repository must be cloned before starting.

---

## Step 1 - Locate your inputs
Confirm you have run Base2Fastq.slurm on the raw AVITI sequencer output, and have access to the following before proceeding:

- The run's `RunManifest.csv`
- The parent `Samples/` directory containing per-sample FASTQ subdirectories

FASTQ files must follow this naming convention:

```
{SampleName}_R1.fastq.gz
{SampleName}_R2.fastq.gz
```

> **Note:** PhiX entries and Unassigned reads are automatically excluded by the pipeline.

---

## Step 2 - Create and activate the conda environment
Run these commands once per installation. Skip if the environment already exists.

```{bash}
# Create the environment from the provided YAML
conda env create -f aviti_read_qc_pipeline.yaml

# Activate it
conda activate aviti_read_qc_pipeline

# Verify a key dependency
snakemake --version   # expected: 9.9.0
```

> **If you already have this environment from an earlier release, update it** — `bwa` and
> `samtools` were added for the optional reference-mapping step:
> ```{bash}
> conda env update -f aviti_read_qc_pipeline.yaml --prune
> ```
> Only needed if you intend to set `mapping.enabled: true`.

---

## Step 3 - Configure the run
Open `config/config.yaml` and fill in the following required fields:

| Parameter | Description |
|---|---|
| `run_name` | Unique identifier for this run (used in output filenames) |
| `run_manifest` | Path to the `RunManifest.csv` |
| `samples_dir` | Path to the parent directory containing FASTQ files |
| `output_dir` | Directory where all outputs will be written |
| `adapter_r1` | Adapter sequence for forward reads |
| `adapter_r2` | Adapter sequence for reverse reads |

Adjust fastp parameters if needed — defaults are appropriate for most AVITI runs.

---

## Step 4 - Dry run (optional: recommended for first time running the pipeline)
Perform a dry run before submitting to SLURM to catch any configuration errors:

```{bash}
conda activate aviti_read_qc_pipeline
cd aviti_read_qc_pipeline
snakemake -n --quiet
```

Check that:
- The printed job list looks correct
- Sample counts match your expectations
- `logs/sample_manifest.log` contains no unexpected grouping warnings

> If sample grouping looks wrong (e.g. unexpected merging of unrelated samples), review the
> Index1 + Index2 pairs in your `RunManifest.csv` before proceeding.

---

## Step 5 - Submit to SLURM

```{bash}
sbatch aviti_read_qc_pipeline.sh
```

Monitor job progress:
```{bash}
squeue -u $USER
```

Per-rule logs are written to `logs/<rule>/`. If a job fails, check the relevant log file there first.

> **Before submitting:** ensure the `source` line inside `aviti_read_qc_pipeline.sh` correctly
> points to your `conda.sh`, and that your NHM email address is set for SLURM notifications.

---

## Step 6 - Review outputs
Once the pipeline completes, review the following:

**SLURM run log:**
```
aviti_read_qc_pipeline/aviti_read_qc_pipeline_{job_id}.err
```
> Does the 'rule all' show 100% completion?

**MultiQC report** — open in a browser:
```
multiqc_report/{run_name}_multiqc_report.html
```
> Does the file exist? Does it open without any errors?

**Per-sample fastp summary CSV:**
```
02_fastp/{run_name}_fastp_summary.csv
```
> Does the file exist? Is it populated with data correctly?

---

## Step 7 - Archive
Once the run is finished (in approximetly 1-5 hours) and inputs are confirmed:
```{bash}
# Compress the output directory into a gzipped tarball
tar czf <name_of_directory_to_tar>.tar.gz <name_of_directory_to_tar/>

# Transfer to appropriate project storage
# (update path as needed)
cp <name_of_directory_to_tar>.tar.g /path/to/storage/
```

---

## Quick reference: expected output structure

```
output_dir/
├── 00_lane_merge/          # Lane-concatenated FASTQs
├── 01_pre_qc/              # falco reports (pre-fastp)
├── 02_fastp/               # Trimmed reads, fastp HTML/JSON, summary CSV
├── 03_post_qc/             # falco reports (post-fastp)
├── 04_seqkit/              # seqkit stats per sample
├── 05_mapping/             # flagstat per sample + summary TSV (only if mapping.enabled)
├── multiqc_report/         # Aggregated MultiQC HTML report
└── logs/                   # Per-rule logs and sample_manifest.log
```

---

## Troubleshooting

| Problem | Likely cause | Fix |
|---|---|---|
| FASTQ files not found | Naming convention mismatch | Check files are named `{Sample}_R1.fastq.gz` |
| Unexpected sample grouping | Duplicate Index1+Index2 pairs | Review `RunManifest.csv`; check `sample_manifest.log` |
| conda activation fails in SLURM | Wrong `source` path | Edit the `source` line in `aviti_read_qc_pipeline.sh` |
| Dry run shows 0 jobs | Config path errors | Check all paths in `config/config.yaml` exist |
| MultiQC report missing samples | Rule failed upstream | Check `logs/` for the failing rule |
| `mapping.reference not found` on startup | `mapping.enabled: true` with a bad or unset `reference` | Correct the path in `config/config.yaml`, or set `mapping.enabled: false` |
| `bwa index` fails writing to the reference directory | Reference is on read-only or shared storage | Build the index yourself where you have write access: `bwa index /path/to/reference.fasta` |
| No Samtools section in the MultiQC report | Mapping disabled, or all flagstat files empty | Confirm `mapping.enabled: true` and check `logs/reference_mapping/` |
| Mapping jobs hit the SLURM walltime | Mapping is hours per sample | Use a longer partition and/or raise `rules.mapping.threads` |

---
