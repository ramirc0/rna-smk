# rna-smk

A [Snakemake](https://github.com/snakemake/snakemake) pipeline for paired-end
RNA-seq.

## Description

Reads are trimmed with fastp and aligned with STAR. The alignments feed bigWig
signal tracks, RSEM gene and isoform quantification and strandedness
inference. A QC report grades every sample on ENCODE's portal audits. Four pixi
workspaces split the launcher from the production, QC and report rule envs, and
profiles cover local and SLURM execution.

### Layout

```
config/
  config.yaml.template   parameters + paths; copy to <run_id>/config.yaml
  samples.tsv.template   one row per sample: sample_id, r1, r2
  <run_id>/              one run's config.yaml and samples.tsv (gitignored)
workflow/
  Snakefile              includes the rule files, defines rule all
  rules/
    common.smk           config, sample sheet, path constants
    preprocess.smk       fastp
    align.smk            STAR index and alignment
    signals.smk          bigWig signal tracks
    gene_quantification.smk  RSEM index and quantification
    qc.smk               gene model, strandedness
    report.smk           gene table, QC table, figures, QC report
  report/qc.yaml         QC report spec: metrics and ENCODE audits
  scripts/               argparse scripts, each with build_parser() + main()
  envs/pixi.toml         production rule env (pixi workspace)
  envs/qc/pixi.toml      QC rule env (a separate pixi workspace)
  envs/report/pixi.toml  QC report env (a separate pixi workspace)
.tests/integration/      ENCODE end-to-end test: config + sample sheet
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env + lint task
pyproject.toml           ruff + numpydoc config
```

### Steps

Preprocessing and alignment, in `preprocess.smk` and `align.smk`:

- `fastp` trims adapters and low-quality bases. It writes HTML and JSON reports.
- `index` builds the STAR genome index with the GTF's splice junctions and
  `sequencing.sjdb_overhang`.
- `align` runs STAR with the ENCODE parameters. It writes a coordinate-sorted
  genome BAM and an RSEM-ready transcriptome BAM (`_anno.bam`), with flagstat
  and STAR log reports as text and JSON.

Signal tracks, in `signals.smk`:

- `bam_to_signals` runs STAR's `inputAlignmentsFromBAM` mode on the genome BAM.
  It writes bigWigs of unique and of all reads, one per strand for stranded
  libraries.

Quantification, in `gene_quantification.smk`:

- `rsem_index` builds the RSEM reference from the FASTA and GTF.
- `rsem_quant` runs `rsem-calculate-expression` with credibility intervals. It
  counts the genes above 1 TPM.

QC, in `qc.smk`:

- `gene_model` converts the GTF to a BED12 gene model with the UCSC tools.
- `strandedness` runs RSeQC's `infer_experiment.py` on the genome BAM. A
  binomial test on the strand fractions gives `forward`, `reverse` or
  `unstranded`. The result is a report only. `sequencing.strand_direction`
  still drives RSEM and the signal tracks. It also keeps the three strand
  fractions in a TSV.

The QC report, in `report.smk`, uses [qc-report](https://github.com/ramirc0/qc-report):

- `gene_table` lists each GTF gene's type and chromosome.
- `gather_qc` collects every sample's QC into one table: fastp, STAR, flagstat,
  RSEM counts and fragment lengths, the share of expected counts by gene type
  (mitochondrial, rRNA, protein-coding, lncRNA) and the strand fractions.
- `qc_figures` draws each sample's read quality, mapping, fragment lengths, gene
  types, expression and strandedness.
- `qc_report` grades each sample on ENCODE's bulk RNA-seq portal audits, read
  depth and read length, and flags cohort outliers. The spec is
  `workflow/report/qc.yaml`.

### Outputs

Each sample writes to `<outdir>/<run_id>/<sample>/<step>/`. The STAR index,
RSEM index, gene model and gene table do not depend on the run. They go to
`<outdir>/refs/`. The QC report goes to `<outdir>/<run_id>/report/`. Logs and benchmarks mirror the run layout under
`logs/<run_id>/` and `benchmarks/<run_id>/`.

| Step directory | Key outputs |
|---|---|
| `processed/` | Trimmed FASTQs, HTML and JSON reports |
| `align/` | Genome BAM, transcriptome BAM, flagstat TXT and JSON, STAR `Log.final.out` and JSON |
| `signals/` | bigWigs: `plusAll`, `plusUniq`, `minusAll`, `minusUniq`, or `all` and `uniq` |
| `rsem_quant/` | `.genes.results`, `.isoforms.results`, genes-detected JSON |
| `strandedness/` | `forward`, `reverse` or `unstranded` in a TXT, strand fractions TSV |
| `report/` (per run) | `qc_report.html`, `qc_metrics.tsv`, `samples/<sample>.js` figures |

### Known issues

These legacy bugs stay to keep outputs identical to the production results.
A fix changes outputs.

- `get_strandedness.py` runs its binomial test on millions of reads, so almost
  any imbalance counts as stranded.

### Conventions

The workflow follows the
[Nextstrain Snakemake style guide](https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html)
and [Snakemake's best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html):

- No `run:` blocks. Custom logic lives in `workflow/scripts/`, called from `shell:`.
- Raw, triple-quoted `shell` blocks (`r"""`), one option per line.
- Every interpolation uses `:q` quoting.
- Config reaches `shell` through `params:` (a `lambda`), never directly.
- `input:` and `output:` are literal path strings. All paths are relative.
- Rules tee logs to file and terminal (`exec &> >(tee {log:q})`).
- Shell blocks end with `sync {output:q}`. On a shared filesystem, Snakemake's
  first stat of an output from the submit host can otherwise wait for the
  writing node to flush it.
- Scripts are called as `{workflow.basedir:q}/scripts/x.py`, so the workflow
  runs from any `--directory`.
- Every rule has `log:` and `benchmark:`. Each names its pixi env in `software:`.
- Docstrings follow numpydoc. `pixi run lint` checks the scripts with ruff and
  `numpydoc lint`. Neither tool parses `.smk` files.

## Getting Started

### Dependencies

- Linux x86-64, the only platform in the pixi workspaces.
- [pixi](https://pixi.sh) on `PATH`.
- Network access to GitHub. Snakemake and its pixi and SLURM plugins install
  from git, and RSEM builds from source.
- A SLURM cluster for `profiles/slurm`.
- A GRCh38 FASTA and a GENCODE GTF. See Inputs.

### Installing

```bash
git clone https://github.com/ramirc0/rna-smk.git
cd rna-smk
pixi install                                                # launcher env
mkdir -p config/<run_id>
cp config/config.yaml.template config/<run_id>/config.yaml  # then edit
cp config/samples.tsv.template config/<run_id>/samples.tsv  # then edit
```

The launcher env (`pixi.toml` at the root) has Snakemake main and its plugins
for pixi deployment and SLURM. The rule envs are three more pixi workspaces.
`workflow/envs/` runs the production rules. RSEM is the one exception to
bioconda. Every bioconda build of RSEM compiles out its `assert()` calls, and
some of them read input, so the credibility intervals go wrong.
`workflow/envs/rsem/` builds RSEM from source with pixi-build instead.
`workflow/envs/qc/` runs the QC rules and `workflow/envs/report/` the QC report.
Each `pixi.lock` pins every package. An `onstart` hook installs every env on
the submit host before the first job. Parallel SLURM jobs would
otherwise race to install them, and the RSEM build clones GitHub, which some
compute nodes cannot reach.

Add a dependency with `pixi add <pkg>` in the workspace's directory. After a
hand edit, run `pixi lock` there (rules use `locked=True`).

### Inputs

The sample sheet is a TSV with a header row. `sample_id` names the sample. `r1`
and `r2` are its raw paired-end FASTQs. Lines that start with `#` are comments.
Each `sample_id` MUST be unique.

The references live in `resources/` (gitignored). Neither may be gzipped.

- the GRCh38 FASTA, no-alt analysis set
- a GENCODE GTF, e.g. `gencode.v49.basic.annotation.gtf`

### Configuration

| Key | Meaning |
|---|---|
| `samples` | Sample sheet path. |
| `outdir` | Results root. |
| `run_id` | Run name. Outputs go under `<outdir>/<run_id>/`. |
| `pixi_env` | Environment of `workflow/envs/pixi.toml` that the production rules use. |
| `references.fasta` | Reference genome FASTA. |
| `references.gtf` | Gene annotation GTF, for STAR, RSEM and the strandedness gene model. |
| `sequencing.is_stranded` | Write plus and minus strand signal tracks. |
| `sequencing.strand_direction` | `forward`, `reverse` or `unstranded`. Sets RSEM's `--forward-prob`. |
| `sequencing.sjdb_overhang` | STAR's `--sjdbOverhang`: read length minus 1. `<outdir>/refs/` holds one index, so a run with another value MUST use its own `outdir`. |
| `sequencing.rsem_seed` | RSEM random seed. |

### Executing program

Each run keeps its config and sample sheet in `config/<run_id>/`. Name the run
`<start date>-<slug>`, e.g. `2026-11-13-encode-k562`. Set `run_id` and
`samples` in the config to match.

```bash
pixi shell                                                  # launcher env
snakemake -n -p --profile profiles/local --configfile config/<run_id>/config.yaml  # dry run
snakemake --profile profiles/local --configfile config/<run_id>/config.yaml        # local
snakemake --profile profiles/slurm --configfile config/<run_id>/config.yaml        # SLURM
pixi run lint                                               # ruff + numpydoc
```

Snakemake MUST run from the repository root. The SLURM profile holds Weng Lab
partitions. Adapt them for another cluster with the
[SLURM executor docs](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
Tune `set-resources` from the `benchmarks/<run_id>/` TSVs after a first run.

`rsem_quant` runs with 22 threads. RSEM's Gibbs sampling and credibility
intervals change with `-p`, so every sample of a cohort MUST use the same
count. The expected counts, TPM and FPKM do not change with `-p`. With fewer
than 22 cores, Snakemake scales the threads down and changes the Gibbs output.

### Testing

`.tests/integration/` runs four ENCODE total RNA-seq replicates end to end on
full GRCh38. They are two replicates each of GM12878 (ENCSR000AEC) and K562
(ENCSR000AEL). The inputs live on the Weng Lab cluster, so it runs only there.

```bash
snakemake --profile profiles/slurm --configfile .tests/integration/config.yaml
```

Outputs go under `results/integration/`. The test builds its own STAR and RSEM
indexes there, so it never touches `results/refs/`. Compare the QC with ENCODE's
for the same replicates.

## Help

- Every rule reran. Any change to a workspace's manifest or lock reruns every
  rule that uses it. So does moving the workdir, because the env hash includes
  the workspace's absolute path.
- A job failed. The profiles set `--show-failed-logs`. The full log sits under
  `logs/<run_id>/`.
- A SLURM node keeps failing jobs. List it under `slurm-exclude-failed-nodes`
  in `profiles/slurm/config.yaml`.

If you have any questions or would like to provide constructive feedback, please open an issue or reach out to the MOHD DACC.

## Authors

Christian Ramirez ([@ramirc0](https://github.com/ramirc0))

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).

## Acknowledgments

Chen, Shifu. 2025. “Fastp 1.0: An Ultra‐fast All‐round Tool for FASTQ Data Quality Control and Preprocessing.” iMeta 4 (5). https://doi.org/10.1002/imt2.70078.

Danecek, Petr, James K. Bonfield, Jennifer Liddle, et al. 2021. “Twelve Years of SAMtools and BCFtools.” GigaScience 10 (2). https://doi.org/10.1093/gigascience/giab008.

Dobin, Alexander, Carrie A. Davis, Felix Schlesinger, et al. 2013. “STAR: Ultrafast Universal RNA-Seq Aligner.” Bioinformatics (Oxford, England) 29 (1): 15–21. https://doi.org/10.1093/bioinformatics/bts635.

ENCODE DCC. 2022. “ENCODE RNA-seq Pipeline,” version 1.2.4. GitHub. https://github.com/ENCODE-DCC/rna-seq-pipeline.

Li, Bo, and Colin N. Dewey. 2011. “RSEM: Accurate Transcript Quantification from RNA-Seq Data with or without a Reference Genome.” BMC Bioinformatics 12 (1). https://doi.org/10.1186/1471-2105-12-323.

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. 2021. “Sustainable Data Analysis with Snakemake.” F1000Research 10 (April): 33. https://doi.org/10.12688/f1000research.29032.2.

Wang, Liguo, Shengqin Wang, and Wei Li. 2012. “RSeQC: Quality Control of RNA-Seq Experiments.” Bioinformatics (Oxford, England) 28 (16): 2184–85. https://doi.org/10.1093/bioinformatics/bts356.
