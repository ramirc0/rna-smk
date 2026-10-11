# CLAUDE.md

RNA-seq Snakemake workflow on the snakemake-template layout.
README.md holds the layout, steps, outputs, and style conventions. Follow it.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi run snakemake -n -p --profile profiles/local
pixi run snakemake --profile profiles/local
pixi run snakemake --profile profiles/slurm --configfile .tests/integration/config.yaml  # ENCODE end-to-end test
```

## Three pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, plus the `lint` env.
- `workflow/envs/pixi.toml` is the production rule env (`SOFTWARE_ENV`).
- `workflow/envs/qc/pixi.toml` is the QC rule env (`QC_ENV`). QC dependencies go there.
- `pyproject.toml` is ruff and numpydoc config only. pixi ignores it.
- Rules use `locked=True`. After a manifest edit, run `pixi lock` in that workspace.
- Any manifest or lock change reruns every rule that uses the workspace. The plugin hashes the whole `pixi.toml` and `pixi.lock`. A QC dependency in `workflow/envs/` would rerun every alignment and quantification.
- `workflow/envs/rsem/` is a pixi-build package: RSEM `8bc1e21` from source. The lock does not hash its recipe. Install it with `pixi install` in `workflow/envs/` before a run in a fresh checkout. Compute nodes may not resolve github.com. After a recipe edit, delete `workflow/envs/.pixi` and rerun the RSEM rules with `--forcerun`.
- Both env paths resolve relative to the rule file. New rule files MUST sit in `workflow/rules/`.

## Outputs must not drift

- `workflow/envs/` pins match `clarity001/rna-smk:latest` as cached on 2025-07-10, the image behind the production results. bedtools cannot match: it is 2.31.1, because 2.30.0 needs zlib <1.3.
- Never use bioconda RSEM. Its builds define `NDEBUG`, which drops the `assert(fscanf(...))` calls that read the Gibbs count vectors. The credibility intervals come out wrong and the step runs for hours.
- Tag `pre-template` is the last commit before the template port. It is the equivalence baseline.
- README.md "Known issues" lists legacy bugs kept on purpose. A fix changes outputs and needs the user's approval.

## Gotchas

- RSEM output depends on `-p`. Keep `rsem_quant` at 44 threads, and `--cores` at 44 or more.
- `<outdir>/refs/` holds one STAR index, built with `sequencing.sjdb_overhang`. A run with another value MUST use its own `outdir`, as the integration test does.
- STAR and RSEM need `ulimit -n 65536`.
- `SIGNAL_TRACKS` MUST NOT use `multiext()`. `expand()` rejects multiext-flagged paths.
- `.smk` files are not linted. Keep helper docstrings numpydoc by hand.
- Resources MUST be flat values. Never derive `mem_mb` or `runtime` from input size.
- The `* attempt` retry scaling in `profiles/slurm/config.yaml` is intentional.
- CLI `--set-threads` replaces the profile's whole mapping. Pass every rule.
- Production lives in `/zata/zippy/ramirezc/rna-smk` and `/data/zusers/ramirezc/rna-smk`. Treat both as read-only. The legacy config's `workdir` points there, and a legacy rerun deletes its STAR and RSEM indexes.

## Verifying a change

No test suite. Diff against the legacy pipeline on real samples:

1. Take the first 1M read pairs of two production samples. Use the full GRCh38 FASTA and GENCODE GTF. STAR's default `--genomeSAindexNbases` is too large for a chr21+chr22 genome.
2. Run `git archive pre-template` with `chmod +x` on its scripts, Snakemake 9.16.3 and the cached image as `--apptainer-prefix`. Give it a scratch `config.yml` whose `workdir` is the scratch dir.
3. Run this workflow with the same per-rule threads.
4. Run `fastp` and `align` single-threaded in both. Multithreading reorders reads.
5. Compare BAM records with `samtools view`, gzipped files decompressed, and bigWigs decoded with `bigtools bigwigtobedgraph`.
