# AVITI_read_QC_pipeline
A snakemake pipeline for QC of raw, basecalled and demultiplexed AVITI24 sequence data, written for the MBL team @NHMUK.

The pipeline takes an AVITI24 RunManifest.csv and a parent directory of raw FASTQ files (Samples/), merges replicate samples across lanes, concatenates replicates where needed, runs pre-QC FastQC (falco), fastp QC, post-QC FastQC (falco) and Seqkit stats, and then produces a summary spreadsheet of fastp QC metrics and aggregates everything into a single MultiQC report. PhiX entries and Unassigned reads are excluded.

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
- fastp=1.3.1
- falco=1.2.5
- multiqc=1.33
- python=3.12
- snakemake=9.9.0
- snakemake-executor-plugin-slurm=1.6.1
- zip=3.0
- seqkit=2.13.0
3. You are now ready to prepare the pipeline to run (see below).

---

## Quick start
1. Follow the installation and conda env creation steps above.
2. Populate `config/config.yaml` with the required run paremeters and paths.
3. Edit conda 'source' line to correctly point to your `conda.sh` and your NHM email address within `aviti_read_qc_pipeline.sh`
4. Run `sbatch aviti_read_qc_pipeline.slurm` to execute the pipeline. This script will submit all required jobs to a SLURM HPC.

> **A detailed, step-by-step SOP can be found [here](https://github.com/SchistoDan/AVITI_read_QC_pipeline/blob/main/SOP_README.md).**

---

## Workflow overview
**The pipeline comprises the following main steps:**
1. The necessary information within the RunManifest.csv is parsed (extracts the [SETTINGS] and [SAMPLES] blocks, while skipping PhiX entries).
2. Samples are then grouped by Index1+Index2 pair to identify replicates across lanes (if any), and per-lane filepaths are constructed for each sample.
3. Together, these form the sample_table dictionary used for downstream rules, and which is visualised in the output sample_manifest.log
4. Rule 1: `lane_merge` - Concatenation of lane replicates (if required). Single-lane samples are copied directly. Empty inputs produce placeholder files that propagate gracefully through downstream rules (calls `workflow/scripts/lane_merge.py`).
5. Rule 2: `pre_fastqc` -  Runs falco implementation of FastQC for R1 and R2 files on the concatenated reads, writing HTML/data/summary files to `01_pre_qc/{sample}/` and zipping them for `multiqc`.
6. Rule 3: `fastp` -  Trims user-specified adapters, filters by quality/length, and deduplicates, writing trimmed reads to `02_fastp/{sample}/` and generating a HTML and JSON report for `multiqc` and `fastp_summary`. Poly-G and poly-X tail trimming, and read error correction, are optional, and along with additional fastp arguments, can be specified in the `config/config.yaml` 
7. Rule 4: `post_fastqc` - Repeats the falco QC but on the fastp-trimmed reads, writing to `03_post_qc/{sample}/` and zipping them for `multiqc`.
8. Rule 5: `seqkit_stats` - Runs `seqkit stats --all --tabular` on both trimmed R1 and R2, writing a tab-separated stats file to `04_seqkit/{sample}/{sample}_seqkit_stats.txt` for `multiqc`.
9. Rule 6: `fastp_summary` - Walks `02_fastp/`, finding all per-sample JSON reports and compiles them into a single CSV at `02_fastp/{run_name}_fastp_summary.csv`. Columns include pre/post read counts, Q20/Q30 rates, GC content, duplication rate, insert size peak, and filtering outcome counts (calls `workflow/scripts/parse_fastp_stats.py`).
10. Rule 6: `multiqc` - Searches `01_pre_qc/`, `02_fastp/`, `03_post_qc/`, and `04_seqkit/` and aggregates all falco zips, fastp JSONs, and seqkit stats files into a single HTML report in `multiqc_report/`.


<div align="center">
  <img width="384" height="546" src="https://github.com/user-attachments/assets/42288bdb-5083-407f-bec8-d51e8da3867a">
</div>



### Sample grouping and lane concatenation
Samples are grouped by matching Index1 + Index2 pairs listed in RunManifest.csv. Samples sharing the same index pair (i.e. the same library sequenced across multiple lanes) are concatenated before QC. A name-based validation layer checks that grouped sample names share a common prefix, and a warning is written to sample_manifest.log if this looks suspicious (but the pipeline does not fail). 

The base sample name used for output files is derived from the longest common prefix (LCP) of the grouped sample names, with trailing non-alphanumeric characters stripped - Examples:
| Grouped names       | Base name |
|---------------------|-----------|
| Pan1, Pan1a         | Pan1      |
| Pan10, Pan10a       | Pan10     |
| AMB04_A07 | AMB04_A07 |
| BGE_0001_A1_1, BGE_0001_A1_2, BGE_0001_A1_3 | BGE_0001_A1 |



### Key parameters in the config.yaml
The following parameters are configurable within the `config/config.yaml`
- `run_name`: Unique identifier for this run (used in output filenames and MultiQC report). Can be anything sensible
- `run_manifest`: Path to the RunManifest.csv file
- `samples_dir`: Path to parent directory under which FASTQ files will be recursively searched. The pipeline expects files named '{SampleName}_R1.fastq.gz' and '{SampleName}_R2.fastq.gz' in per-sample subdirectories anywhere under this root
- `output_dir`: Directory where all pipeline outputs will be written (created if absent)
- `adapter_r1` & `adapter_r2`: Adapter sequences used for forward and reverse reads
- Fastp:
  - `qualified_quality_phred`: Phred score threshold; bases below this are 'unqualified'
  - `unqualified_percent_limit`: Maximum % of unqualified bases allowed per read before it is discarded
  - `dedup`: Flag to set whether identical read pairs are deduplicated, or not
  - `trim_poly_g`: If enabled, Fastp can detect the polyG in read tails and trims them. PolyG can happen in read tails since G means no signal in the Illumina two-color systems
  - `trim_poly_x`: If enabled, Fastp can detect polyX (e.g. polyA) tails and trims them. If polyG and polyX tail trimming are both enabled, fastp will perform polyG trimming first, then perform polyX trimming
  - `correction`: Overlap-based base correction for paired-end reads. it has the following parameters by default: overlap_len_require (default 30), overlap_diff_limit (default 5) and overlap_diff_percent_limit (default 20%)
  - `extra_args`: Any additional fastp arguments that are not already configurable
- MultiQC:
  - `extra_args`: Any additional multiqc arguments that are not already configurable
- rules: Controls memory and threading resource allocation given to each rule. 


---

## Output directory structure
```
output_dir/
├── 00_lane_merge/          # Per-sample lane-concatenated FASTQ files
│   └── {sample}/
│       ├── {sample}_R1.fastq.gz
│       └── {sample}_R2.fastq.gz
├── 01_pre_qc/              # FastQC reports on merged (pre-fastp) reads
│   └── {sample}/
│       ├── {sample}_R1_fastqc.html
│       ├── {sample}_R1_fastqc.zip
│       ├── {sample}_R2_fastqc.html
│       └── {sample}_R2_fastqc.zip
├── 02_fastp/               # fastp-trimmed reads and QC reports
│   └── {sample}/
│       ├── {sample}_R1.fastq.gz
│       ├── {sample}_R2.fastq.gz
│       ├── {sample}_fastp.html
│       └── {sample}_fastp.json
├── 03_post_qc/             # FastQC reports on fastp-trimmed reads
│   └── {sample}/
│       ├── {sample}_R1_fastqc.html
│       ├── {sample}_R1_fastqc.zip
│       ├── {sample}_R2_fastqc.html
│       └── {sample}_R2_fastqc.zip
├── multiqc_report/
│   ├── {run_name}_multiqc_report.html
│   └── {run_name}_multiqc_data/
└── logs/
    ├── sample_manifest.log     # Which samples were processed, grouped, and skipped
    ├── lane_merge/
    ├── pre_qc/
    ├── fastp/
    ├── post_qc/
    └── multiqc/
```

---

## Benchmarking
End-to-end, the pipeline ran on 96 low coverage WGS (i.e.genome skims) generated from museum specimens, run on both flowcell lanes (i.e. replicates across lanes), in **4 hours and 56 minutes** with modest resources < lane_merge = 16GB/8 threads; fastqc = 8GB/4 threads; fastp = 16GB/8 threads; seqkit = 8GB/4 threads; multiqc = 16GB/2 threads.

---

## Citations & authorship
This snakemake pipeline was written by Dan Parsons for NHMUK Molecualr Biology Laboratories.

Please see below for a list of citations for tools utilised in the pipeline:
| Tool | URL | citation | Version |
| --- | --- | --- | --- |
| snakemake | https://snakemake.readthedocs.io/en/stable/ | [Mölder et al., 2025](https://f1000research.com/articles/10-33/v3) | 9.9.0 | 
| Falco | https://github.com/smithlabcode/falco | [DS Brandine and Smith, 2021](https://f1000research.com/articles/8-1874/v2) | 1.2.5 |
| Fastp | https://github.com/opengene/fastp | [Chen, 2025](https://onlinelibrary.wiley.com/doi/10.1002/imt2.70078) | 1.3.1 |
| MultiQC | https://github.com/MultiQC/MultiQC | [Ewels et al., 2016](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507) | 1.33 |
| Seqkit | https://github.com/shenwei356/seqkit | [Shen et al., 2024](https://onlinelibrary.wiley.com/doi/10.1002/imt2.191) | 2.13.0 |

