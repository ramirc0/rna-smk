"""Write each gene's ID, name, type and chromosome from a GENCODE GTF."""

import argparse

import polars as pl


def build_parser():
    """Build the gene_table.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the gene_table.py flags.
    """
    p = argparse.ArgumentParser(description="Tabulate a GTF's genes.")
    p.add_argument("--gtf", required=True, help="GENCODE GTF, not gzipped.")
    p.add_argument("--output", required=True, help="TSV to write.")
    return p


def main(argv=None):
    """Write the gene table: ``gene_id``, ``gene_name``, ``gene_type`` and ``chrom``.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    attributes = pl.col("attributes")
    (
        pl.scan_csv(
            args.gtf,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            quote_char=None,
            new_columns=["chrom", "source", "feature", "start", "end"]
            + ["score", "strand", "frame", "attributes"],
            schema_overrides={"chrom": pl.String, "attributes": pl.String},
        )
        .filter(pl.col("feature") == "gene")
        .select(
            gene_id=attributes.str.extract(r'gene_id "([^"]+)"'),
            gene_name=attributes.str.extract(r'gene_name "([^"]+)"'),
            gene_type=attributes.str.extract(r'gene_type "([^"]+)"'),
            chrom="chrom",
        )
        .sink_csv(args.output, separator="\t")
    )


if __name__ == "__main__":
    main()
