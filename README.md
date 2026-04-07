# AVITI_read_QC_pipeline
A snakemake pipeline for QC of raw, baecalled and demultiplexed AVITI24 sequence data, compilaiton and visualisation of pre- and post-QC metrics, and concatenation of duplicate sequencing attempts across flowcell lanes

Takes a the AVITI24 RunManifest.csv and a parent directory of FASTQ files (Samples/), merges samples across lanes, concatenates where needed, runs pre-QC FastQC, fastp QC, post-QC FastQC, and aggregates everything into a single MultiQC report. PhiX entries and Unassigned reads are excluded.



## Dependencies & installation
1. Clone this repository.
2. Install all necessary dependencies listed in `aviti_read_qc_pipeline.yaml` provided in this repository, using the following commands:
```
# conda must be installed first
conda env create -f aviti_read_qc_pipeline.yaml

# Check the environment was succesfully created
conda activate aviti_read_qc_pipeline
```

**Core dependencies:**
- fastp=1.3.1
- falco=1.2.5
- multiqc=1.33
- python=3.12
- snakemake=9.9.0
- snakemake-executor-plugin-slurm=1.6.1
3. You are now ready to prepare the pipeline to run (see below).



## Quick start
1. Follow the installation and conda env creation steps above.
2. Populate `config/config.yaml` with the required run paremeters and paths.
3. Run `aviti_read_qc_pipeline.slurm` with `sbatch` to execute the pipeline. This script will submit all required jobs to a SLURM HPC.



## Workflow overview
**Rule worklow:** lane_merge → pre_fastqc → fastp → post_fastqc → multiqc


### Sample grouping and lane concatenation
Samples are grouped by matching Index1 + Index2 pairs listed in RunManifest.csv. Samples sharing the same index pair (i.e. the same library sequenced across multiple lanes) are concatenated before QC. A name-based validation layer checks that grouped sample names share a common prefix, and a warning is written to sample_manifest.log if this looks suspicious (but the pipeline does not fail). 

The base sample name used for output files is derived from the longest common prefix (LCP) of the grouped sample names, with trailing non-alphanumeric characters stripped - Examples:
| Grouped names       | Base name |
|---------------------|-----------|
| Pan1, Pan1a         | Pan1      |
| Pan10, Pan10a       | Pan10     |
| AMB04_A07 | AMB04_A07 |
| BGE_0001_A1_1, BGE_0001_A1_2, BGE_0001_A1_3 | BGE_0001_A1 |



### Key parameter explanations
- `qualified_quality_phred`: defines what "unqualified" means in terms of base quality (phred) score
- `unqualified_percent_limit`: defines what the 'tolerance' is 
- `trim_poly_g`: 
- `trim_poly_x`: 

### Python helper functions
- validate_config(config) — Sanity-checks the config.yaml before anything else runs. Confirms all required keys exist, that path values point to real files/directories, and that numeric fastp parameters are valid types. Exits with a clear error message rather than letting Snakemake fail cryptically mid-run.
parse_run_manifest(manifest_path) — Reads the RunManifest.csv line by line, handling the two-section format ([SETTINGS] and [SAMPLES]), skipping comment lines and blank lines. Returns the settings dict, the list of sample rows, and a list of already-skipped entries (PhiX).
group_samples_by_index(samples_rows) — Takes the flat list of sample rows and groups them into a dict keyed by (Index1, Index2) tuple. This is the primary grouping authority — samples sharing an index pair are the same biological library across lanes.
longest_common_prefix(strings) — Pure string utility. Compares a list of strings character by character and returns the longest prefix shared by all of them. Used to derive base sample names from grouped lane names.
derive_base_name(names) — Calls longest_common_prefix then strips any trailing non-alphanumeric characters (e.g. trailing underscores or hyphens) from the result. Converts ["Pan1", "Pan1a"] → "Pan1" cleanly.
validate_name_grouping(base_name, names, index_key) — The name-matching validation layer on top of index grouping. Warns (but never fails) if the derived base name is suspiciously short (< 3 chars) or if any sample name in the group doesn't start with the base name. Results go to the manifest log.
find_fastq_files(sample_name, samples_dir) — Recursively globs under samples_dir for files matching {sample_name}_*_R1*.fastq.gz and _R2*. Applies strict_filter to prevent prefix collisions (e.g. Pan1 matching Pan1a files) by anchoring the regex to require an underscore immediately after the sample name.
build_sample_table(config) — Orchestrates everything above. Calls parse → group → find files for each group → validates names → assembles the final sample_table dict that all Snakemake rules read from. Also accumulates the processed/skipped/warnings lists for logging.
write_manifest_log(...) — Writes a human-readable summary of the entire sample resolution process to logs/sample_manifest.log before any jobs run. Critically useful for verifying grouping is correct before committing cluster resources.


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

## Citations & authorship
This snakemake pipeline was written by Dan Parsons for NHMUK Molecualr Biology Laboratories.

Please see below for a list of citations for tools utilised in the pipeline:
| Tool | URL | citation | Version |
| --- | --- | --- | --- |
| snakemake | ... | 9.9.0 | 
| Falco | https://github.com/smithlabcode/falco | [de Sena Brandine and Smith, 2021](https://f1000research.com/articles/8-1874/v2) | 1.2.5 |
| Fastp | https://github.com/opengene/fastp | [Chen, 2025](https://onlinelibrary.wiley.com/doi/10.1002/imt2.70078) | 1.3.1 |
| MultiQC | ... | ... | 1.33 |
