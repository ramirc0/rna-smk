"""Collect every sample's QC into one table, a row per sample, for the run's report."""

import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import polars as pl
from qc_outputs import (
    composition,
    read_fragment_lengths,
    read_gene_counts,
    read_rsem_counts,
    read_star_log,
)
from qc_report import read_fastp


def build_parser():
    """Build the gather_qc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the gather_qc.py flags.
    """
    p = argparse.ArgumentParser(description="Collect each sample's QC metrics.")
    p.add_argument("--run-dir", required=True, help="Run directory of sample dirs.")
    p.add_argument("--samples", required=True, nargs="+", help="Sample IDs.")
    p.add_argument("--gene-table", required=True, help="TSV from gene_table.py.")
    p.add_argument("--threads", type=int, default=1, help="Samples read at once.")
    p.add_argument("--output", required=True, help="TSV to write.")
    return p


def sample_metrics(run, sample, gene_table):
    """Collect one sample's metrics, keyed by the ids in ``workflow/report/qc.yaml``.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    sample : str
        Sample ID.
    gene_table : polars.DataFrame
        Gene types and chromosomes, from gene_table.py.

    Returns
    -------
    dict
        One row of the metrics table.
    """
    d = run / sample
    fastp_json = d / f"processed/{sample}.json"
    fastp = read_fastp(fastp_json)
    summary = json.loads(fastp_json.read_text())["summary"]
    star = read_star_log(d / f"align/{sample}_Log.final.json")
    flagstat = json.loads((d / f"align/{sample}_genome_flagstat.json").read_text())
    rsem = read_rsem_counts(d / f"rsem_quant/{sample}.stat/{sample}.cnt")
    lengths = read_fragment_lengths(d / f"rsem_quant/{sample}.stat/{sample}.model")
    genes_detected = json.loads(
        (d / f"rsem_quant/{sample}_no_genes_detected.json").read_text()
    )
    strand = pl.read_csv(d / f"strandedness/{sample}_fractions.tsv", separator="\t")
    genes = read_gene_counts(d / f"rsem_quant/{sample}.genes.results", gene_table)
    return {
        "sample_id": sample,
        "read_pairs": fastp["reads"] // 2,
        "read_length": fastp["read_length"],
        "q30_pct": fastp["q30_pct"],
        "gc_pct": fastp["gc_pct"],
        "duplication_pct": fastp["duplication_pct"],
        "passed_filter_pct": 100
        * summary["after_filtering"]["total_reads"]
        / summary["before_filtering"]["total_reads"],
        "read_depth": star["Uniquely mapped reads number"]
        + star["Number of reads mapped to multiple loci"],
        "unique_pct": star["Uniquely mapped reads %"],
        "multi_pct": star["% of reads mapped to multiple loci"],
        "too_many_loci_pct": star["% of reads mapped to too many loci"],
        "too_short_pct": star["% of reads unmapped: too short"],
        "mismatch_rate_pct": star["Mismatch rate per base, %"],
        "annotated_splices_pct": 100
        * star["Number of splices: Annotated (sjdb)"]
        / star["Number of splices: Total"],
        "mapped_length": star["Average mapped length"],
        "properly_paired_pct": flagstat["samtools_genome_flagstat"][
            "paired_properly_pct"
        ],
        # RSEM counts fragments; STAR counts read pairs as one input read.
        "transcriptome_pct": 100 * rsem["alignable"] / star["Number of input reads"],
        "genes_detected": genes_detected["number_of_genes_gt_1_TPM"],
        "mean_fragment_length": (lengths["insert_size"] * lengths["probability"]).sum(),
        **composition(genes),
        "strand_forward_pct": 100 * strand["forward"][0],
        "strand_reverse_pct": 100 * strand["reverse"][0],
        "strand_undetermined_pct": 100 * strand["undetermined"][0],
    }


def main(argv=None):
    """Write the run's metrics table.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    run = Path(args.run_dir)
    gene_table = pl.read_csv(
        args.gene_table, separator="\t", columns=["gene_id", "gene_type", "chrom"]
    )
    with ThreadPoolExecutor(args.threads) as pool:
        rows = list(
            pool.map(lambda s: sample_metrics(run, s, gene_table), args.samples)
        )
    pl.DataFrame(rows).write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
