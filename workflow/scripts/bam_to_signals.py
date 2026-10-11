"""Run the bam to signals (bigWig) step of the ENCODE rna-seq-pipeline."""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import logging
import os
import subprocess

logger = logging.getLogger(__name__)

TRACKS = {
    "stranded": {
        "Signal.UniqueMultiple.str1.out.bg": "minusAll",
        "Signal.Unique.str1.out.bg": "minusUniq",
        "Signal.UniqueMultiple.str2.out.bg": "plusAll",
        "Signal.Unique.str2.out.bg": "plusUniq",
    },
    "unstranded": {
        "Signal.UniqueMultiple.str1.out.bg": "all",
        "Signal.Unique.str1.out.bg": "uniq",
    },
}


def run(command):
    """Log and run a command, raising on a nonzero exit.

    Parameters
    ----------
    command : list of str
        Program and arguments.
    """
    logger.info("Running %s", " ".join(command))
    subprocess.run(command, check=True)


def build_parser():
    """Build the bam_to_signals.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the bam_to_signals.py flags.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", required=True, help="Output directory.")
    parser.add_argument(
        "--threads", type=int, default=1, help="Threads for bigWig conversion."
    )
    parser.add_argument("--bamfile", required=True, help="Coordinate-sorted BAM.")
    parser.add_argument(
        "--chrom_sizes", required=True, help="Chromosome sizes of the BAM's genome."
    )
    parser.add_argument(
        "--strandedness",
        required=True,
        choices=list(TRACKS),
        help="`stranded` writes one track per strand.",
    )
    parser.add_argument(
        "--bamroot",
        required=True,
        help="Output name root. `foo` writes foo_plusAll.bw and so on.",
    )
    return parser


def main(argv=None):
    """Write the STAR signal tracks of a BAM as bigWigs.

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
    out = args.output_dir
    run(
        [
            "STAR",
            "--runMode", "inputAlignmentsFromBAM",
            "--inputBAMfile", args.bamfile,
            "--outWigType", "bedGraph",
            "--outWigStrand", args.strandedness.capitalize(),
            "--outWigReferencesPrefix", "chr",
            "--outFileNamePrefix", f"{out}/",
        ]
    )  # fmt: skip
    for bedgraph, track in TRACKS[args.strandedness].items():
        bedgraph = os.path.join(out, bedgraph)
        # STAR writes chromosomes in genome order, not the lexical order that
        # bigtools' default `--sorted all` expects.
        run(
            [
                "bigtools", "bedgraphtobigwig",
                "--sorted", "start",
                "--nthreads", str(args.threads),
                bedgraph,
                args.chrom_sizes,
                os.path.join(out, f"{args.bamroot}_{track}.bw"),
            ]
        )  # fmt: skip
        os.remove(bedgraph)


if __name__ == "__main__":
    main()
