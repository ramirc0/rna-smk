"""Read the QC files this pipeline writes, for the run's report."""

import json
from pathlib import Path

import polars as pl

# GENCODE gene types whose RNA is ribosomal.
RRNA_TYPES = ("rRNA", "rRNA_pseudogene", "Mt_rRNA")


def read_star_log(path):
    """Read a STAR ``Log.final.json`` into numbers.

    Parameters
    ----------
    path : str or pathlib.Path
        ``<sample>_Log.final.json``, written by align.py.

    Returns
    -------
    dict of str to float
        Each STAR metric with numbers, keyed as STAR names it. Percent signs
        are dropped, so ``"91.15%"`` reads as 91.15. Dates are left out.
    """
    log = json.loads(Path(path).read_text())["star_log_qc"]
    numbers = {}
    for name, value in log.items():
        try:
            numbers[name] = float(value.rstrip("%"))
        except ValueError:
            continue
    return numbers


def read_rsem_counts(path):
    """Read the read counts from RSEM's ``.cnt`` file.

    Parameters
    ----------
    path : str or pathlib.Path
        ``<sample>.stat/<sample>.cnt``.

    Returns
    -------
    dict of str to int
        ``unalignable``, ``alignable``, ``filtered`` and ``total`` reads,
        from the file's first line.
    """
    first = Path(path).read_text().split("\n", 1)[0].split()
    unalignable, alignable, filtered, total = map(int, first)
    return {
        "unalignable": unalignable,
        "alignable": alignable,
        "filtered": filtered,
        "total": total,
    }


def read_fragment_lengths(path):
    """Read RSEM's fragment length distribution from its ``.model`` file.

    Parameters
    ----------
    path : str or pathlib.Path
        ``<sample>.stat/<sample>.model`` of a paired-end run.

    Returns
    -------
    polars.DataFrame
        ``insert_size`` and ``probability``. RSEM writes ``lb ub span`` on the
        fifth line and the probability of each length in ``(lb, ub]`` on the
        sixth (PairedEndModel::write, LenDist::write).
    """
    lines = Path(path).read_text().splitlines()
    low, high, _ = map(int, lines[4].split())
    probabilities = [float(p) for p in lines[5].split()]
    return pl.DataFrame(
        {"insert_size": range(low + 1, high + 1), "probability": probabilities}
    )


def read_gene_counts(genes_results, gene_table):
    """Join RSEM's gene counts and TPM with each gene's type and chromosome.

    Parameters
    ----------
    genes_results : str or pathlib.Path
        RSEM's ``<sample>.genes.results``.
    gene_table : polars.DataFrame
        ``gene_id``, ``gene_type`` and ``chrom`` per gene, from gene_table.py.

    Returns
    -------
    polars.DataFrame
        ``gene_id``, ``expected_count``, ``TPM``, ``gene_type`` and ``chrom``.
    """
    counts = pl.read_csv(
        genes_results,
        separator="\t",
        columns=["gene_id", "expected_count", "TPM"],
    )
    return counts.join(gene_table, on="gene_id", how="left")


def composition(genes):
    """Share of RSEM's expected counts in each class of gene, in percent.

    Parameters
    ----------
    genes : polars.DataFrame
        Output of `read_gene_counts`.

    Returns
    -------
    dict of str to float
        ``mito_pct`` (genes on chrM), ``rrna_pct`` (`RRNA_TYPES`),
        ``protein_coding_pct`` and ``lncrna_pct``.
    """
    count = pl.col("expected_count")
    return (
        genes.select(
            mito_pct=count.filter(pl.col("chrom") == "chrM").sum(),
            rrna_pct=count.filter(pl.col("gene_type").is_in(RRNA_TYPES)).sum(),
            protein_coding_pct=count.filter(
                pl.col("gene_type") == "protein_coding"
            ).sum(),
            lncrna_pct=count.filter(pl.col("gene_type") == "lncRNA").sum(),
        )
        .select(pl.all() * 100 / genes["expected_count"].sum())
        .row(0, named=True)
    )
