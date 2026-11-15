"""Write each sample's report figures, from the data each QC step keeps."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import polars as pl
from qc_outputs import (
    read_fragment_lengths,
    read_gene_counts,
    read_rsem_counts,
    read_star_log,
)
from qc_report import insert_size_figure, read_quality_figure, write_sample_panels

# STAR's fate of each read pair, in the order the bar stacks them.
STAR_FATES = {
    "Unique": "Uniquely mapped reads %",
    "Multi": "% of reads mapped to multiple loci",
    "Too many loci": "% of reads mapped to too many loci",
    "Too short": "% of reads unmapped: too short",
    "Mismatches": "% of reads unmapped: too many mismatches",
    "Other": "% of reads unmapped: other",
    "Chimeric": "% of chimeric reads",
}
# Gene types drawn on their own. The rest pool as "other".
GENE_TYPES = 8


def build_parser():
    """Build the qc_figures.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the qc_figures.py flags.
    """
    p = argparse.ArgumentParser(description="Write each sample's report figures.")
    p.add_argument("--run-dir", required=True, help="Run directory of sample dirs.")
    p.add_argument("--report-dir", required=True, help="Directory of the report.")
    p.add_argument("--samples", required=True, nargs="+", help="Sample IDs.")
    p.add_argument("--gene-table", required=True, help="TSV from gene_table.py.")
    p.add_argument("--threads", type=int, default=1, help="Samples drawn at once.")
    return p


def mapping_figure(star):
    """Draw where STAR put the read pairs, as one stacked bar.

    Parameters
    ----------
    star : dict of str to float
        Output of `qc_outputs.read_star_log`.

    Returns
    -------
    plotly.graph_objects.Figure
        One segment per fate in `STAR_FATES`, in percent of input pairs.
    """
    fig = go.Figure(
        [
            go.Bar(
                x=[star[key]],
                y=["Read pairs"],
                name=fate,
                orientation="h",
                hovertemplate=f"{fate}: %{{x:.2f}}%<extra></extra>",
            )
            for fate, key in STAR_FATES.items()
        ]
    )
    fig.update_layout(barmode="stack", xaxis_title="Read pairs (%)", height=220)
    fig.update_yaxes(showticklabels=False)
    return fig


def gene_type_figure(genes):
    """Draw the share of RSEM's expected counts per gene type.

    Parameters
    ----------
    genes : polars.DataFrame
        Output of `qc_outputs.read_gene_counts`.

    Returns
    -------
    plotly.graph_objects.Figure
        The `GENE_TYPES` largest types, largest first, and the rest as "other".
    """
    shares = (
        genes.group_by("gene_type")
        .agg(pl.col("expected_count").sum())
        .with_columns(
            share=100 * pl.col("expected_count") / pl.col("expected_count").sum()
        )
        .sort("share", descending=True)
    )
    top = shares.head(GENE_TYPES)
    types = top["gene_type"].to_list() + ["other"]
    values = top["share"].to_list() + [shares["share"].tail(-GENE_TYPES).sum()]
    fig = go.Figure(
        go.Bar(
            x=values,
            y=types,
            orientation="h",
            hovertemplate="%{y}: %{x:.2f}%<extra></extra>",
        )
    )
    fig.update_yaxes(autorange="reversed")
    fig.update_layout(xaxis_title="Expected counts (%)", showlegend=False)
    return fig


def expression_figure(genes):
    """Draw the distribution of gene TPM, on a log scale.

    Parameters
    ----------
    genes : polars.DataFrame
        Output of `qc_outputs.read_gene_counts`.

    Returns
    -------
    plotly.graph_objects.Figure
        Histogram of log10 TPM over expressed genes, with 1 TPM marked.
    """
    # Binned here, so the page holds 60 bars instead of every gene's TPM.
    tpm = genes.filter(pl.col("TPM") > 0)["TPM"].to_numpy()
    counts, edges = np.histogram(np.log10(tpm), bins=60)
    fig = go.Figure(
        go.Bar(
            x=(edges[:-1] + edges[1:]) / 2,
            y=counts,
            width=np.diff(edges),
            hovertemplate="log10 TPM %{x:.2f}<br>%{y:,} genes<extra></extra>",
        )
    )
    fig.add_vline(x=0, line_dash="dot", annotation_text="1 TPM")
    fig.update_layout(xaxis_title="log10 TPM", yaxis_title="Genes", showlegend=False)
    return fig


def strand_figure(fractions):
    """Draw infer_experiment.py's strand fractions.

    Parameters
    ----------
    fractions : polars.DataFrame
        The strandedness rule's ``_fractions.tsv``: ``forward``, ``reverse``
        and ``undetermined``.

    Returns
    -------
    plotly.graph_objects.Figure
        One bar per fraction, in percent of sampled reads.
    """
    fig = go.Figure(
        go.Bar(
            x=fractions.columns,
            y=[100 * fractions[c][0] for c in fractions.columns],
            hovertemplate="%{x}: %{y:.1f}%<extra></extra>",
        )
    )
    fig.update_layout(yaxis_title="Reads (%)", yaxis_range=[0, 100], showlegend=False)
    return fig


def write_panels(sample, run, report_dir, gene_table):
    """Write one sample's ``samples/<sample>.js``.

    Parameters
    ----------
    sample : str
        Sample ID.
    run : pathlib.Path
        Run directory.
    report_dir : pathlib.Path
        Directory of the report HTML.
    gene_table : polars.DataFrame
        Gene types and chromosomes, from gene_table.py.
    """
    d = run / sample
    stat = d / f"rsem_quant/{sample}.stat"
    lengths = read_fragment_lengths(stat / f"{sample}.model")
    alignable = read_rsem_counts(stat / f"{sample}.cnt")["alignable"]
    genes = read_gene_counts(d / f"rsem_quant/{sample}.genes.results", gene_table)
    write_sample_panels(
        report_dir / "samples" / f"{sample}.js",
        sample,
        [
            ("Read quality", read_quality_figure(d / f"processed/{sample}.json")),
            (
                "Mapping",
                mapping_figure(read_star_log(d / f"align/{sample}_Log.final.json")),
            ),
            (
                "Fragment lengths",
                insert_size_figure(
                    lengths.select(
                        "insert_size",
                        fragments=(pl.col("probability") * alignable).round(),
                    )
                ),
            ),
            ("Gene types", gene_type_figure(genes)),
            ("Expression", expression_figure(genes)),
            (
                "Strandedness",
                strand_figure(
                    pl.read_csv(
                        d / f"strandedness/{sample}_fractions.tsv", separator="\t"
                    )
                ),
            ),
        ],
    )


def main(argv=None):
    """Write every sample's figures.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    gene_table = pl.read_csv(
        args.gene_table, separator="\t", columns=["gene_id", "gene_type", "chrom"]
    )
    draw = partial(
        write_panels,
        run=Path(args.run_dir),
        report_dir=Path(args.report_dir),
        gene_table=gene_table,
    )
    with ProcessPoolExecutor(args.threads) as pool:
        list(pool.map(draw, args.samples))


if __name__ == "__main__":
    main()
