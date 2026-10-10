"""Infer library strandedness from a BAM with RSeQC's `infer_experiment.py`."""

import argparse
import logging
import re
import subprocess
import sys

from scipy.stats import binomtest


def build_parser():
    """Build the get_strandedness.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the get_strandedness.py flags.
    """
    parser = argparse.ArgumentParser(
        description="Parse the output of RSeQC's `infer_experiment.py` to determine "
        "strandedness. Writes 'forward', 'reverse', or 'unstranded'.",
    )
    parser.add_argument("--bam", required=True, help="Path to input bam.")
    parser.add_argument(
        "--gene_model", required=True, help="BED12 gene model for infer_experiment.py."
    )
    parser.add_argument("--output", required=True, help="Text file to write.")
    return parser


def main(argv=None):
    """Write `forward`, `reverse` or `unstranded` for a BAM.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    logger = logging.getLogger(__name__)

    logger.info("Counting reads in input bam...")
    samtools_count = subprocess.run(
        ["samtools", "view", "-c", "-q", "31", args.bam],
        capture_output=True,
        text=True,
        check=True,
    )
    total_reads = int(samtools_count.stdout.strip())

    logger.info("Running RSeQC infer_experiment.py...")
    infer_experiment = subprocess.run(
        ["infer_experiment.py", "-i", args.bam, "-r", args.gene_model],
        capture_output=True,
        text=True,
        check=True,
    )
    lines = [line.strip() for line in infer_experiment.stdout.splitlines()[2:]]
    endedness_pattern = r"This is (.*) Data"
    fraction_pattern = r"Fraction of reads failed to determine: (.*)"
    forward_pattern = r'Fraction of reads explained by "1\+\+,1--,2\+-,2-\+": (.*)'
    reverse_pattern = r'Fraction of reads explained by "1\+-,1-\+,2\+\+,2--": (.*)'

    logger.info("Parsing RSeQC infer_experiment.py output...")
    forward = None
    reverse = None
    try:
        for line in lines:
            if re.match(endedness_pattern, line):
                _ = re.search(endedness_pattern, line).group(1)
            elif re.match(fraction_pattern, line):
                _ = float(re.search(fraction_pattern, line).group(1))
            elif re.match(forward_pattern, line):
                forward = float(re.search(forward_pattern, line).group(1))
            elif re.match(reverse_pattern, line):
                reverse = float(re.search(reverse_pattern, line).group(1))
    except Exception:
        logger.exception("Could not parse RSeQC infer_experiment.py output.")
        sys.exit(1)

    forward_reads = int(total_reads * forward)
    reverse_reads = int(total_reads * reverse)
    determinable_reads = forward_reads + reverse_reads

    try:
        p_value = binomtest(forward_reads, determinable_reads, p=0.5).pvalue
    except Exception:
        logger.exception("Could not calculate p-value.")
        sys.exit(1)

    logger.info(f"P-value: {p_value:f}")
    if p_value < 0.05:
        logger.info("Library is stranded")
        strandedness = "forward" if forward_reads > reverse_reads else "reverse"
    else:
        logger.info("Library is unstranded")
        strandedness = "unstranded"
    with open(args.output, "w") as f:
        print(strandedness, file=f)


if __name__ == "__main__":
    main()
