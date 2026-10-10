"""Infer library strandedness from a BAM with RSeQC's `infer_experiment.py`."""

import argparse
import gzip
import logging
import os
import re
import subprocess
import sys
import tempfile

import requests
from scipy.stats import binomtest


def build_parser():
    """Build the get_strandedness.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the get_strandedness.py flags.
    """
    parser = argparse.ArgumentParser(
        prog="get-strandedness",
        description="Takes a BAM file as input and parses the output of RSeQC's `infer_experiment.py` to determine strandedness. Returns 'forward', 'reverse', or 'unstranded' to stdout.",
    )
    parser.add_argument("--src", help="View source code", action="store_true")
    parser.add_argument(
        "input",
        help="Path to input bam. Pass '-' to read from stdin.",
        nargs="?",
        metavar="<IN.BAM>",
        type=str,
    )
    return parser


def main(argv=None):
    """Print `forward`, `reverse` or `unstranded` for a BAM.

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

    if args.src:
        try:
            editor = os.environ.get("EDITOR")
            if not editor:
                editor = "vi"
            subprocess.check_call([editor, sys.argv[0]])
            sys.exit(0)
        except subprocess.CalledProcessError:
            print("Error: Could not open file in editor.")
    in_bam = args.input

    with (
        tempfile.NamedTemporaryFile(delete=False) as temp1,
        tempfile.NamedTemporaryFile(delete=False) as temp2,
    ):
        logger.info("Counting reads in input bam...")
        samtools_count = subprocess.run(
            ["samtools", "view", "-c", "-q", "31", in_bam],
            capture_output=True,
            text=True,
            check=True,
        )
        total_reads = int(samtools_count.stdout.strip())

        logger.info(
            "Downloading latest gene model from 'https://sourceforge.net/projects/rseqc/files/latest/download'..."
        )
        response = requests.get(
            "https://sourceforge.net/projects/rseqc/files/latest/download",
            allow_redirects=True,
        )
        temp1.write(gzip.decompress(response.content))

        latest_genemodel = temp1.name

        logger.info("Running RSeQC infer_experiment.py...")
        subprocess.check_call(
            ["infer_experiment.py", "-i", in_bam, "-r", latest_genemodel], stdout=temp2
        )

        with open(temp2.name) as f:
            lines = f.readlines()

        lines = [line.strip() for line in lines[2:]]
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
            if forward_reads > reverse_reads:
                print("forward")
            else:
                print("reverse")
        else:
            logger.info("Library is unstranded")
            print("unstranded")

    os.unlink(temp1.name)
    os.unlink(temp2.name)


if __name__ == "__main__":
    main()
