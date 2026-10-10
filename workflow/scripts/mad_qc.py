"""Compare the gene FPKMs of two RSEM runs by the MAD of their log ratios.

A port of the madqc Rust tool. It keeps its math, including the upper middle
value as the median of an even count and the uncentered standard deviation.
"""

import argparse
import json

import matplotlib.pyplot as plt
import numpy as np
import poikilos as pk
import polars as pl


def build_parser():
    """Build the mad_qc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the mad_qc.py flags.
    """
    p = argparse.ArgumentParser(
        description="Calculate the MAD of log ratios of two RSEM genes.results files."
    )
    p.add_argument("--tsv1", required=True, help="RSEM genes.results of one sample.")
    p.add_argument("--tsv2", required=True, help="RSEM genes.results of the other.")
    p.add_argument("--json", required=True, help="Metrics JSON to write.")
    p.add_argument(
        "--plot", required=True, help="MA plot SVG to write (PNG beside it)."
    )
    return p


def mad_qc(fpkm1, fpkm2):
    """Compute the MAD QC metrics and the MA values of two FPKM vectors.

    Parameters
    ----------
    fpkm1 : numpy.ndarray
        FPKM per gene of one sample, with no zeros.
    fpkm2 : numpy.ndarray
        FPKM of the same genes in the other sample, with no zeros.

    Returns
    -------
    metrics : dict of str to float
        `MAD`, `Standard_deviation` and `Pearson_correlation` of the genes
        with a mean log2 FPKM above 0.
    a : numpy.ndarray
        Mean log2 FPKM per gene.
    m : numpy.ndarray
        Log2 FPKM ratio per gene.
    """
    log1 = np.log2(fpkm1)
    log2 = np.log2(fpkm2)
    a = (log1 + log2) / 2.0
    m = log1 - log2
    mask = a > 0

    m_kept = m[mask]
    mad = float(np.sort(np.abs(m_kept))[m_kept.size // 2]) * 1.4826
    sd = float(np.sqrt(np.mean(m_kept**2)))
    corr = float(np.corrcoef(log1[mask], log2[mask])[0, 1])
    metrics = {"MAD": mad, "Standard_deviation": sd, "Pearson_correlation": corr}
    return metrics, a, m


def main(argv=None):
    """Join two genes.results files, write the metrics JSON and the MA plot.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    df1 = pl.read_csv(args.tsv1, separator="\t", infer_schema_length=None)
    df2 = pl.read_csv(args.tsv2, separator="\t", infer_schema_length=None)
    df = (
        df1.join(df2, on="gene_id", how="inner", maintain_order="left")
        .select("gene_id", "TPM", "FPKM", "TPM_right", "FPKM_right")
        .filter(~((pl.col("FPKM") == 0.0) | (pl.col("FPKM_right") == 0.0)))
    )

    metrics, a, m = mad_qc(df["FPKM"].to_numpy(), df["FPKM_right"].to_numpy())
    print(f"MAD of log ratios: {metrics['MAD']}")
    print(f"Standard deviation of log ratios: {metrics['Standard_deviation']}")
    print(f"Pearson correlation: {metrics['Pearson_correlation']}")
    with open(args.json, "w") as f:
        json.dump(metrics, f, indent=2)

    pk.use("plain")
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(a, m, s=2, rasterized=True)
    ax.set(xlim=(-10, 20), ylim=(-10, 10), xlabel="A", ylabel="M", title="MAD QC Plot")
    pk.despine(ax)
    pk.save_figure(fig, args.plot)
    plt.close(fig)


if __name__ == "__main__":
    main()
