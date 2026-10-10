# RNA-seq Snakemake Pipeline

## Overview

Here we present a [snakemake](https://github.com/snakemake/snakemake.git) pipeline for RNA-seq preprocessing, alignment, signal generation, gene quantification, and quality control across multiple samples.

> [!NOTE]
> Alignment and quantification steps are adapted from the [ENCODE RNA-seq pipeline](https://github.com/ENCODE-DCC/rna-seq-pipeline). See individual script headers for attribution details.

## Example Usage

### Prerequisites

This pipeline expects paired-end (or single-end) FASTQ files as input. You will need a metadata TSV file describing your samples and their corresponding read files (see [Example Configfile](#example-configfile) below).

### Dependencies

This pipeline requires the following dependencies:

- [`conda` (preferably `mamba`)](https://github.com/conda-forge/miniforge)
- [`singularity`](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)
- (Optional) [`slurm`](https://slurm.schedmd.com/quickstart.html)

To start, first setup the runtime environment:

```bash
mamba env create --name=rna-smk --file=environments/rna-smk.yaml
mamba activate rna-smk
```

### Example Configfile

For a complete example, see [config/config.yml](config/config.yml). From this file you'll want to pay attention to the following sections:

```yaml
workdir: "/path/to/rna-smk"
```

- `workdir`: The absolute path to the working directory for the pipeline. All relative paths in the config are resolved from here.

```yaml
sequencing:
  endedness: "paired"
  is_stranded: True
  strand_direction: "reverse"
  rsem_seed: 42
```

- `sequencing.endedness`: Whether the library is `"paired"` or `"single"` ended.
- `sequencing.is_stranded`: Whether the library is stranded (`True` or `False`). This controls whether strand-specific signal tracks (bigWig files) are generated.
- `sequencing.strand_direction`: The strandedness of the library. Must be one of `"forward"`, `"reverse"`, or `"unstranded"`. If you are unsure, the included utility script `workflow/scripts/get_strandedness.py` can infer strandedness from a BAM file.
- `sequencing.rsem_seed`: Random seed for RSEM quantification reproducibility.

```yaml
input_files:
  metadata: "jobs/2026-01-06/metadata.tsv"
  gtf_file: "jobs/2026-01-06/gencode.v49.basic.annotation.gtf"
  reference_fasta: "jobs/2026-01-06/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
```

- `input_files.metadata`: Path to the sample metadata file. This is a tab-delimited file with a header containing at least three columns: `Sample`, `R1`, and `R2`. Lines beginning with `#` are treated as comments. `Sample` values *must be unique*. `R1` and `R2` contain absolute paths to the corresponding FASTQ files.
- `input_files.gtf_file`: Path to the gene annotation GTF file. Must **not** be gzipped.
- `input_files.reference_fasta`: Path to the reference genome FASTA file. Must **not** be gzipped.

```yaml
container: docker://clarity001/rna-smk:latest
```

- `container`: The Docker/Singularity container image used for most pipeline rules. This image bundles STAR, RSEM, samtools, fastp, and bigtools.

### Defining Resources

All resource definitions for rules are located in [profile/slurm/config.yaml](profile/slurm/config.yaml).

```yaml
  align:
    slurm_extra: "--cpu-freq=High-High:Performance"
    threads: 44
    cpus_per_task: 44
    mem_gb: 250
    mem_mb: 250000
```

In this instance, we are requesting a job on a SLURM cluster with 44 threads, 44 CPUs per task, and 250 GB of memory.

These resource definitions are specific to the Weng Lab's SLURM cluster. If you are running the pipeline on your own cluster, you'll need to adjust these values accordingly. We recommend consulting the [snakemake-executor-plugin-slurm documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) for further information.

### Running the Pipeline

With the hard part out of the way, you can now run the pipeline:

```bash
snakemake --workflow-profile profile/slurm
```

> [!NOTE]
> It is recommended that you run the pipeline in dry run mode first (add the `-n` flag). Also note that `snakemake` must be ran in the root of this repository.

This will run the pipeline and produce all outputs to a directory `results`.

### Pipeline Steps

The pipeline executes the following steps in order:

1. **Preprocessing** (`fastp`): Adapter trimming and quality filtering of raw FASTQ files.
2. **Indexing** (`index`): Builds a STAR genome index, with splice junction annotation from the provided GTF.
3. **Alignment** (`align`): Aligns reads to the genome and generates both a genome-coordinate BAM and a transcriptome-coordinate (annotation) BAM. When using STAR, alignment retries up to 3 times with progressively longer SLURM partitions.
4. **Signal Generation** (`bam_to_signals`): Converts genome BAMs to bigWig signal tracks. Produces strand-specific tracks (plus/minus, all/unique) for stranded libraries, or unstranded tracks otherwise.
5. **RSEM Indexing** (`rsem_index`): Builds an RSEM reference from the genome FASTA and GTF.
6. **RSEM Quantification** (`rsem_quant`): Quantifies gene and isoform expression from the annotation BAM using RSEM.

### Outputs

The pipeline produces results organized under the `results/` directory:

| Directory | Key Outputs |
|---|---|
| `results/processed/` | Trimmed FASTQs, HTML and JSON QC reports from fastp |
| `results/index/` | STAR genome index, chromosome sizes file |
| `results/align/{sample}/` | Genome-coordinate BAM, transcriptome-coordinate (annotation) BAM, flagstat JSON, STAR `Log.final.out` |
| `results/signals/{sample}/` | Strand-specific bigWig signal tracks (plus/minus, all/unique) or unstranded bigWig tracks |
| `results/rsem_index/` | RSEM reference index |
| `results/rsem_quant/{sample}` | Gene and isoform quantification tables (`.genes.results`, `.isoforms.results`), gene detection JSON |

### What if I cannot run on a SLURM cluster?

If you cannot run the pipeline on a SLURM cluster, you can run it locally. However, this is not a supported use case.

To run locally, go to the [workflow/Snakefile](workflow/Snakefile) and add a `localrules` directive listing the rules you wish to run on the local machine. Note that alignment and indexing steps are memory-intensive and may not be feasible on a local workstation.

### Utility Scripts

The `workflow/scripts/` directory contains a utility script:

- `get_strandedness.py`: Infers library strandedness from a BAM file using RSeQC's `infer_experiment.py`. Returns `"forward"`, `"reverse"`, or `"unstranded"` to stdout.

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

> Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. _STAR: ultrafast universal RNA-seq aligner_. Bioinformatics, 29(1), 15–21, **2013**. https://doi.org/10.1093/bioinformatics/bts635.

> Li, B. & Dewey, C. N. _RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome_. BMC Bioinformatics, 12, 323, **2011**. https://doi.org/10.1186/1471-2105-12-323.

## Questions

If you have any questions or would like to provide constructive feedback, please open an issue or reach out to the MOHD DACC.
