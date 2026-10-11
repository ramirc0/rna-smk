"""Run the RSEM quantification step of the ENCODE rna-seq-pipeline."""

import argparse
import json
import logging
import subprocess

import pandas as pd

logger = logging.getLogger(__name__)

STRAND_TO_FORWARD_PROB = {"forward": 1, "unstranded": 0.5, "reverse": 0}


def build_parser():
    """Build the rsem.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the rsem.py flags.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--rsem_index", type=str, help="path to rsem index")
    parser.add_argument("--anno_bam", type=str, help="STAR alignment to annotation.")
    parser.add_argument(
        "--read_strand", type=str, choices=["forward", "reverse", "unstranded"]
    )
    parser.add_argument("--rnd_seed", type=int, help="random seed", default=42)
    parser.add_argument("--ncpus", type=int, help="number of cpus available")
    parser.add_argument("--mem_mb", type=int, help="memory available in GB")
    parser.add_argument("--outprefix", type=str, help="prefix for output files")
    return parser


def main(argv=None):
    """Quantify with RSEM and count the genes above 1 TPM.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s: %(message)s",
    )
    rsem_call = [
        "rsem-calculate-expression",
        "--bam",
        "--estimate-rspd",
        "--calc-ci",
        "--seed", str(args.rnd_seed),
        "-p", str(args.ncpus),
        "--no-bam-output",
        "--ci-memory", str(args.mem_mb),
        "--forward-prob", str(STRAND_TO_FORWARD_PROB[args.read_strand]),
        "--paired-end",
        args.anno_bam,
        args.rsem_index,
        args.outprefix,
    ]  # fmt: skip
    logger.info("Running RSEM command %s", " ".join(rsem_call))
    subprocess.check_call(rsem_call)

    quant_tsv = f"{args.outprefix}.genes.results"
    quants = pd.read_csv(quant_tsv, sep="\t", usecols=["TPM"])
    no_genes_detected = sum(quants["TPM"] > 1)
    qc_record = {"number_of_genes_gt_1_TPM": no_genes_detected}
    logger.info("Number of genes greater than 1 TPM: %s", no_genes_detected)
    with open(f"{args.outprefix}_no_genes_detected.json", "w") as f:
        json.dump(qc_record, f, indent=4)


if __name__ == "__main__":
    main()
