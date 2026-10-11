"""Run the RSEM quantification step of the ENCODE rna-seq-pipeline."""

import argparse
import json
import logging
import shlex
import subprocess
from collections import OrderedDict

import pandas as pd

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)

RSEM_COMMAND = """
    rsem-calculate-expression --bam \
    --estimate-rspd \
    --calc-ci \
    --seed {rnd_seed} \
    -p {ncpus} \
    --no-bam-output \
    --ci-memory {mem_mb} \
    --forward-prob {fwd_prob} \
    --paired-end \
    {anno_bam} \
    {resem_index} \
    {outprefix}
    """


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
    strand_to_fwd_prob = {"forward": 1, "unstranded": 0.5, "reverse": 0}
    gene_detection_threshold = 1

    rsem_call = shlex.split(
        RSEM_COMMAND.format(
            rnd_seed=args.rnd_seed,
            ncpus=args.ncpus,
            mem_mb=args.mem_mb,
            fwd_prob=strand_to_fwd_prob[args.read_strand],
            anno_bam=args.anno_bam,
            resem_index=args.rsem_index,
            outprefix=args.outprefix,
        )
    )
    logger.info("Running RSEM command %s", " ".join(rsem_call))
    subprocess.check_call(rsem_call)

    quant_tsv = str(args.outprefix) + ".genes.results"
    quants = pd.read_csv(quant_tsv, sep="\t", usecols=["TPM"])
    no_genes_detected = sum(quants["TPM"] > gene_detection_threshold)
    qc_record = OrderedDict()
    qc_record["number_of_genes_gt_1_TPM"] = no_genes_detected
    logger.info("Number of genes greater than 1 TPM: %s", no_genes_detected)
    with open(str(args.outprefix) + "_no_genes_detected.json", "w") as f:
        json.dump(qc_record, f, indent=4)


if __name__ == "__main__":
    main()
