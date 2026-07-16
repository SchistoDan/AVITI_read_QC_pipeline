<img width="300" height="150" alt="NHM_logo_new" src="https://github.com/user-attachments/assets/7d72e101-621a-4b3b-9d28-2bb7e5e2a085" />

# AVITI_read_QC_pipeline
 
A Snakemake pipeline for QC of raw, basecalled and demultiplexed AVITI24 sequence data, written for the MBL/SeqFac team @NHMUK.
 
The pipeline parses one or more AVITI24 `RunManifest.csv` files and their sibling `Samples/` directories, concatenates lane replicates (where required), runs pre-QC falco (FastQC-compatible), fastp adapter trimming and quality filtering, post-QC falco, and seqkit stats. It then produces a per-run fastp summary CSV and aggregates everything into a single MultiQC report. PhiX entries are excluded automatically.
 
---

## Dependencies & installation
1. Clone this repository.
2. Install all necessary dependencies listed in `aviti_read_qc_pipeline.yaml` provided in this repository, using the following commands:
```
# conda must be installed first
conda env create -f aviti_read_qc_pipeline.yaml

# Check the environment was successfully created
conda activate aviti_read_qc_pipeline
```

**Core dependencies:**
| Tool | Version |
|------|---------|
| python | 3.12 |
| snakemake | 9.9.0 |
| snakemake-executor-plugin-slurm | 1.6.1 |
| fastp | 1.3.1 |
| falco | 1.2.5 |
| multiqc | 1.33 |
| seqkit | 2.13.0 |
| zip | 3.0 |
 
3. You are now ready to configure and run the pipeline (see below).

---

## Quick start
 
1. Follow the installation and conda environment creation steps above.
2. Populate `config/config.yaml` with the required run parameters and paths (see [Key parameters](#key-parameters-in-configconfigyaml) below).
3. Edit the conda `source` line in `aviti_read_qc_pipeline.sh` to point to your `conda.sh`, and update your NHM email address.
4. Run `sbatch aviti_read_qc_pipeline.slurm` to submit all jobs to a SLURM HPC cluster.
> **A detailed, step-by-step SOP can be found [here](https://github.com/NHM-Sequencing-Facility/AVITI_read_QC_pipeline/blob/main/SOP_README.md).**

---

## Workflow overview
The pipeline comprises the following main steps:
1. **Manifest parsing** — One or more `RunManifest.csv` files are parsed, extracting the `[SETTINGS]` and `[SAMPLES]` blocks and skipping PhiX entries. The `Samples/` directory is derived automatically from each manifest's location; no separate `samples_dir` config key is required. FASTQ files are indexed in a single recursive pass rather than per-sample glob calls (important for performance on networked HPC filesystems).
2. **Sample grouping** — Samples are grouped by matching `Index1 + Index2` pair to identify lane replicates. When `additional_run_manifests` are provided, samples appearing in more than one run have their FASTQ file lists merged across runs. Settings blocks from all manifests are compared and the pipeline exits with a clear diff if any key differs. A summary is written to `logs/sample_manifest.log`.
3. **Rule 1 — `lane_merge`** — Lane replicates are concatenated into a single R1/R2 pair per sample (calls `workflow/scripts/lane_merge.py`). Single-lane samples are copied directly. When `lane_merge.enabled: false`, single-lane samples are symlinked instead of copied; multi-lane samples raise an error (merging is required and cannot be skipped). Empty inputs produce placeholder files that propagate gracefully through all downstream rules.
4. **Rule 2 — `pre_fastqc`** — Runs falco on the merged R1 and R2 reads, writing HTML/data/summary files and zip archives to `01_pre_qc/{sample}/` for MultiQC.
5. **Rule 3 — `fastp`** — Trims adapters, filters by quality and length, and deduplicates reads. Poly-G and poly-X tail trimming and overlap-based base correction are optional. Trimmed reads and HTML/JSON reports are written to `02_fastp/{sample}/`.
6. **Rule 4 — `post_fastqc`** — Repeats falco QC on the fastp-trimmed reads, writing to `03_post_qc/{sample}/`.
7. **Rule 5 — `seqkit_stats`** — Runs `seqkit stats --all --tabular` on both trimmed R1 and R2, writing a tab-separated stats file to `04_seqkit/{sample}/{sample}_seqkit_stats.txt`.
8. **Rule 6 — `fastp_summary`** — Walks `02_fastp/`, finds all per-sample JSON reports, and compiles them into a single CSV at `02_fastp/{run_name}_fastp_summary.csv`. Also writes `{run_name}_general_stats_mqc.yaml`, which injects accurate raw and final read/base counts (from `summary.before_filtering` and `summary.after_filtering`) into the MultiQC General Statistics table — correcting for the fastp module's default use of `filtering_result.passed_filter_reads`, which is measured before deduplication (calls `workflow/scripts/parse_fastp_stats.py`).
9. **Rule 7 — `multiqc`** — Searches `01_pre_qc/`, `02_fastp/`, `03_post_qc/`, and `04_seqkit/`, and aggregates all falco zips, fastp JSONs, seqkit stats files, and the custom general-stats YAML into a single HTML report in `multiqc_report/`.


<div align="center">
  <img width="384" height="546" src="https://github.com/user-attachments/assets/42288bdb-5083-407f-bec8-d51e8da3867a">
</div>



### Sample grouping and lane concatenation
Samples are grouped by matching `Index1 + Index2` pairs listed in the `RunManifest.csv`. Samples sharing the same index pair (i.e. the same library sequenced across multiple lanes) are concatenated before QC. A name-based validation layer checks that grouped sample names share a common prefix; a warning is written to `sample_manifest.log` if this looks suspicious, but the pipeline does not fail.
 
The base sample name used for output files is derived from the longest common prefix (LCP) of the grouped sample names, with trailing non-alphanumeric characters stripped. Examples:
 
| Grouped names | Base name |
|---|---|
| Pan1, Pan1a | Pan1 |
| Pan10, Pan10a | Pan10 |
| AMB04_A07 | AMB04_A07 |
| BGE_0001_A1_1, BGE_0001_A1_2, BGE_0001_A1_3 | BGE_0001_A1 |


### Cross-run merging
When `additional_run_manifests` (in the `config.yaml`) is populated, the pipeline merges FASTQ file lists for samples that appear in more than one run. The `[SETTINGS]` blocks from all manifests must be identical; any mismatch causes a hard exit with a diff of the conflicting keys. Samples present in only one run are included as-is. Cross-run merging always requires `lane_merge.enabled: true`, since merged samples will by definition have more than one input file per read direction.


## Key parameters in `config/config.yaml`
**General**
| Parameter | Description |
|---|---|
| `run_name` | Unique identifier for this run, used in output filenames and the MultiQC report title |
| `run_manifest` | Path to the primary `RunManifest.csv`. The sibling `Samples/` directory is derived automatically from this path |
| `additional_run_manifests` | Optional list of additional `RunManifest.csv` paths for cross-run sample merging. Requires `lane_merge.enabled: true` |
| `output_dir` | Directory where all pipeline outputs will be written (created if absent) |
 
**Lane merge** 
| Parameter | Description | Default |
|---|---|---|
| `lane_merge.enabled` | If `true`, concatenate multi-lane samples. If `false`, symlink single-lane samples through without copying; multi-lane samples raise an error. Must be `true` if `additional_run_manifests` is set | `true` |
 
**fastp**
| Parameter | Description | Default |
|---|---|---|
| `adapter_r1` | Adapter sequence for R1 reads | Illumina TruSeq R1 |
| `adapter_r2` | Adapter sequence for R2 reads | Illumina TruSeq R2 |
| `qualified_quality_phred` | Phred score threshold; bases below this are considered unqualified | `15` |
| `unqualified_percent_limit` | Maximum % of unqualified bases per read before the read is discarded | `40` |
| `min_length` | Minimum read length after trimming; shorter reads are discarded | `15` |
| `dedup` | Remove PCR/optical duplicate read pairs | `true` |
| `trim_poly_g` | Trim poly-G tails (recommended for two-colour Illumina/AVITI chemistry) | `true` |
| `trim_poly_x` | Trim poly-X (e.g. poly-A) tails. Applied after poly-G trimming if both are enabled | `false` |
| `correction` | Overlap-based base correction for paired-end reads (overlap_len_require 30, overlap_diff_limit 5, overlap_diff_percent_limit 20%) | `false` |
| `extra_args` | Any additional fastp arguments as a raw string | `""` |
 
**Falco (pre & post QC)**
| Parameter | Description | Default |
|---|---|---|
| `falco.extra_args` | Any additional falco arguments as a raw string (applied to both pre- and post-QC falco runs) | `""` |

**MultiQC**
| Parameter | Description | Default |
|---|---|---|
| `multiqc.extra_args` | Any additional MultiQC arguments as a raw string | `""` |
 
**Resource allocation (`rules`)**
Each rule block accepts `mem_mb`, `threads`, and `partition` (SLURM partition name). Memory is scaled by retry attempt number on failure (up to the configured `retries` count per rule).
| Rule | Default mem_mb | Default threads |
|---|---|---|
| `lane_merge` | 16384 | 8 |
| `fastqc` (pre & post) | 16384 | 8 |
| `fastp` | 16384 | 8 |
| `seqkit` | 16384 | 8 |
| `multiqc` | 16384 | 2 |


---

## Output directory structure
```
output_dir/
├── 00_lane_merge/              # Per-sample lane-concatenated FASTQ files
│   └── {sample}/
│       ├── {sample}_R1.fastq.gz
│       └── {sample}_R2.fastq.gz
├── 01_pre_qc/                  # falco reports on merged (pre-fastp) reads
│   └── {sample}/
│       ├── {sample}_R1_fastqc.html
│       ├── {sample}_R1_fastqc.zip
│       ├── {sample}_R2_fastqc.html
│       └── {sample}_R2_fastqc.zip
├── 02_fastp/                   # fastp-trimmed reads and QC reports
│   ├── {sample}/
│   │   ├── {sample}_R1.fastq.gz
│   │   ├── {sample}_R2.fastq.gz
│   │   ├── {sample}_fastp.html
│   │   └── {sample}_fastp.json
│   ├── {run_name}_fastp_summary.csv      # Compiled per-sample fastp metrics
│   └── {run_name}_general_stats_mqc.yaml # Raw/final read+base counts for MultiQC
├── 03_post_qc/                 # falco reports on fastp-trimmed reads
│   └── {sample}/
│       ├── {sample}_R1_fastqc.html
│       ├── {sample}_R1_fastqc.zip
│       ├── {sample}_R2_fastqc.html
│       └── {sample}_R2_fastqc.zip
├── 04_seqkit/                  # seqkit stats on fastp-trimmed reads
│   └── {sample}/
│       └── {sample}_seqkit_stats.txt
├── multiqc_report/
│   ├── {run_name}_multiqc_report.html
│   └── {run_name}_multiqc_report_data/
└── logs/
    ├── sample_manifest.log     # Samples processed, grouped, skipped, and warnings
    ├── lane_merge/
    ├── pre_qc/
    ├── fastp/
    ├── post_qc/
    ├── seqkit/
    ├── fastp_summary/
    └── multiqc/
```

---

## Benchmarking
End-to-end, the pipeline ran on 96 low-coverage WGS (genome skims) generated from museum specimens, sequenced across both flowcell lanes (i.e. replicates per lane), in **4 hours and 56 minutes** with the following resources: `lane_merge` 16 GB / 8 threads; `fastqc` 8 GB / 4 threads; `fastp` 16 GB / 8 threads; `seqkit` 8 GB / 4 threads; `multiqc` 16 GB / 2 threads.

---

## Citations & authorship
This Snakemake pipeline was written by Dan Parsons for the NHMUK Molecular Biology Laboratories/Sequencing Facility.
 
| Tool | URL | Citation | Version |
|---|---|---|---|
| Snakemake | https://snakemake.readthedocs.io/en/stable/ | [Mölder et al., 2025](https://f1000research.com/articles/10-33/v3) | 9.9.0 |
| Falco | https://github.com/smithlabcode/falco | [Brandine & Smith, 2021](https://f1000research.com/articles/8-1874/v2) | 1.2.5 |
| fastp | https://github.com/OpenGene/fastp | [Chen, 2025](https://onlinelibrary.wiley.com/doi/10.1002/imt2.70078) | 1.3.1 |
| MultiQC | https://github.com/MultiQC/MultiQC | [Ewels et al., 2016](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507) | 1.33 |
| seqkit | https://github.com/shenwei356/seqkit | [Shen et al., 2024](https://onlinelibrary.wiley.com/doi/10.1002/imt2.191) | 2.13.0 |

