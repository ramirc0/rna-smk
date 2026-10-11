# rna-smk

A [Snakemake](https://github.com/snakemake/snakemake) pipeline for paired-end
RNA-seq. Reads are trimmed and aligned with STAR. The alignments feed bigWig
signal tracks, RSEM gene and isoform quantification and strandedness
inference.

> [!NOTE]
> Alignment, signal and quantification steps come from the
> [ENCODE RNA-seq pipeline](https://github.com/ENCODE-DCC/rna-seq-pipeline).
> The script headers give attribution.

## Layout

```
config/
  config.yaml.template   parameters + paths; copy to config.yaml (gitignored)
  samples.example.tsv    one row per sample: sample_id, r1, r2
workflow/
  Snakefile              includes the rule files, defines rule all
  rules/
    common.smk           config, sample sheet, path constants
    preprocess.smk       fastp
    align.smk            STAR index and alignment
    signals.smk          bigWig signal tracks
    gene_quantification.smk  RSEM index and quantification
    qc.smk               gene model, strandedness
  scripts/               argparse scripts, each with build_parser() + main()
  envs/pixi.toml         production rule env (pixi workspace)
  envs/qc/pixi.toml      QC rule env (a separate pixi workspace)
.tests/integration/      ENCODE end-to-end test: config + sample sheet
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env + lint task
pyproject.toml           ruff + numpydoc config
```

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
cp config/samples.example.tsv config/samples.tsv   # then edit
pixi install --manifest-path workflow/envs/pixi.toml  # builds RSEM once
pixi shell                                         # launcher env
snakemake -n -p --profile profiles/local           # dry run
snakemake --profile profiles/local                 # local
snakemake --profile profiles/slurm                 # SLURM
```

Snakemake MUST run from the repository root. The SLURM profile holds Weng Lab
partitions. Adapt them for another cluster with the
[SLURM executor docs](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
Tune `set-resources` from the `benchmarks/<run_id>/` TSVs after a first run.

`rsem_quant` MUST run with 44 threads. RSEM's Gibbs sampling and credibility
intervals change with `-p`, so another thread count gives other numbers. The
local profile sets `cores: 44` for this reason.

## Integration test

`.tests/integration/` runs four ENCODE total RNA-seq replicates end to end on
full GRCh38. They are two replicates each of GM12878 (ENCSR000AEC) and K562
(ENCSR000AEL). The inputs live on the Weng Lab cluster, so it runs only there.

```bash
snakemake --profile profiles/slurm --configfile .tests/integration/config.yaml
```

Outputs go under `results/integration/`. Its reads are 101 bp, so it sets
`sjdb_overhang: 100` and builds its own STAR and RSEM indexes there. It never
touches `results/refs/`. Compare the QC with ENCODE's
for the same replicates.

## Inputs

The sample sheet is a TSV with a header row. `sample_id` names the sample. `r1`
and `r2` are its raw paired-end FASTQs. Lines that start with `#` are comments. Each `sample_id`
MUST be unique.

The references live in `resources/` (gitignored). Neither may be gzipped.

- the GRCh38 FASTA, no-alt analysis set
- a GENCODE GTF, e.g. `gencode.v49.basic.annotation.gtf`

## Configuration

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

## Steps

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
  still drives RSEM and the signal tracks.

## Outputs

Each sample writes to `<outdir>/<run_id>/<sample>/<step>/`. The STAR index,
RSEM index and gene model do not depend on the run. They go to
`<outdir>/refs/`. Logs and benchmarks mirror the run layout under
`logs/<run_id>/` and `benchmarks/<run_id>/`.

| Step directory | Key outputs |
|---|---|
| `processed/` | Trimmed FASTQs, HTML and JSON reports |
| `align/` | Genome BAM, transcriptome BAM, flagstat TXT and JSON, STAR `Log.final.out` and JSON |
| `signals/` | bigWigs: `plusAll`, `plusUniq`, `minusAll`, `minusUniq`, or `all` and `uniq` |
| `rsem_quant/` | `.genes.results`, `.isoforms.results`, genes-detected JSON |
| `strandedness/` | `forward`, `reverse` or `unstranded` in a TXT |

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake main and its plugins
for pixi deployment and SLURM. The rule envs are two more pixi workspaces.
`workflow/envs/` runs the production rules. Its tool versions match
`clarity001/rna-smk:latest` as cached on 2025-07-10, the image that produced the
production results. RSEM is the one exception to bioconda. Every bioconda build
of RSEM compiles out its `assert()` calls, and some of them read input, so the
credibility intervals go wrong. `workflow/envs/rsem/` builds the image's RSEM
commit from source with pixi-build instead. `workflow/envs/qc/` runs the QC rules. Each `pixi.lock` pins
every package. pixi MUST be on `PATH`. Jobs install the envs on first use.
The production env is the exception. Its RSEM build clones GitHub, which some
compute nodes cannot reach, and concurrent jobs would build it at once. Install
it on the login node before the first run.

Add a dependency with `pixi add <pkg>` in the workspace's directory. After a
hand edit, run `pixi lock` there (rules use `locked=True`). Any change to a
workspace's manifest or lock reruns every rule that uses it. So does moving the
workdir, because the env hash includes the workspace's absolute path.

## Conventions

The workflow follows the
[Nextstrain Snakemake style guide](https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html)
and [Snakemake's best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html):

- No `run:` blocks. Custom logic lives in `workflow/scripts/`, called from `shell:`.
- Raw, triple-quoted `shell` blocks (`r"""`), one option per line.
- Every interpolation uses `:q` quoting.
- Config reaches `shell` through `params:` (a `lambda`), never directly.
- `input:` and `output:` are literal path strings. All paths are relative.
- Rules tee logs to file and terminal (`exec &> >(tee {log:q})`).
- Every rule has `log:` and `benchmark:`. Each names its pixi env in `software:`.
- Docstrings follow numpydoc. `pixi run lint` checks the scripts with ruff and
  `numpydoc lint`. Neither tool parses `.smk` files.

## Differences from ENCODE

The STAR alignment flags and the RSEM flags match
[ENCODE's pipeline](https://github.com/ENCODE-DCC/rna-seq-pipeline) v1.2.4.
The steps around them differ:

| Area | ENCODE v1.2.4 | rna-smk |
|---|---|---|
| Trimming | None | fastp |
| STAR, RSEM, samtools | 2.5.1b, 1.2.31, 1.9 | 2.7.11b, 1.3.3, 1.22 |
| STAR index | `--sjdbOverhang 100`, with ERCC and phiX spike-ins | `sequencing.sjdb_overhang`, no spike-ins |
| Annotation | GENCODE V29 primary assembly plus tRNAs | The configured GTF. Production uses GENCODE v49 basic. |
| bigWig conversion | kent `bedSort` and `bedGraphToBigWig` | bigtools |
| RSEM seed | 12345 | `sequencing.rsem_seed`, 42 in the template |
| Genes-detected JSON | `QCMetricRecord` layout | One key, `number_of_genes_gt_1_TPM` |
| Strandedness inference | None | RSeQC, report only |
| Kallisto, reads by gene type, MAD QC, pBAM | Yes | No |
| Several FASTQ pairs per replicate, single-end reads | Yes | No |

The trimming, tool versions and annotation change every output. The bigWig
tools change only the file bytes.

## Known issues

These legacy bugs stay to keep outputs identical to the production results.
A fix changes outputs.

- `get_strandedness.py` runs its binomial test on millions of reads, so almost
  any imbalance counts as stranded.

## References

Chen, Shifu. 2025. “Fastp 1.0: An Ultra‐fast All‐round Tool for FASTQ Data Quality Control and Preprocessing.” iMeta 4 (5). https://doi.org/10.1002/imt2.70078.

Danecek, Petr, James K. Bonfield, Jennifer Liddle, et al. 2021. “Twelve Years of SAMtools and BCFtools.” GigaScience 10 (2). https://doi.org/10.1093/gigascience/giab008.

Dobin, Alexander, Carrie A. Davis, Felix Schlesinger, et al. 2013. “STAR: Ultrafast Universal RNA-Seq Aligner.” Bioinformatics (Oxford, England) 29 (1): 15–21. https://doi.org/10.1093/bioinformatics/bts635.

Li, Bo, and Colin N. Dewey. 2011. “RSEM: Accurate Transcript Quantification from RNA-Seq Data with or without a Reference Genome.” BMC Bioinformatics 12 (1). https://doi.org/10.1186/1471-2105-12-323.

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. 2021. “Sustainable Data Analysis with Snakemake.” F1000Research 10 (April): 33. https://doi.org/10.12688/f1000research.29032.2.

Wang, Liguo, Shengqin Wang, and Wei Li. 2012. “RSeQC: Quality Control of RNA-Seq Experiments.” Bioinformatics (Oxford, England) 28 (16): 2184–85. https://doi.org/10.1093/bioinformatics/bts356.

## Questions

If you have any questions or would like to provide constructive feedback, please open an issue or reach out to the MOHD DACC.
