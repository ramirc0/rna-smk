"""Run the bam to signals (bigWig) step of the ENCODE rna-seq-pipeline."""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import logging
import os
import shlex
import subprocess
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)

STAR_COMMAND = """STAR --runMode inputAlignmentsFromBAM \
                --inputBAMfile {input_bam} \
                --outWigType bedGraph \
                --outWigStrand {strandedness} \
                --outWigReferencesPrefix chr
                --outFileNamePrefix {output_dir}/"""


def build_parser():
    """Build the bam_to_signals.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the bam_to_signals.py flags.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_dir", type=str, help="Path to output directory.", required=True
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="Number of threads for `bigtools bedgraphtobigwig`",
        default=1,
    )
    parser.add_argument("--bamfile", type=str, help="Input bam")
    parser.add_argument(
        "--chrom_sizes", type=str, help="chromosome sizes file the input bam"
    )
    parser.add_argument("--is_stranded", type=bool)
    parser.add_argument(
        "--bamroot",
        type=str,
        help="""
             Root name for output bams. For example out_bam
             will create out_bam_genome.bam and out_bam_anno.bam
             """,
        default="out_bam",
    )
    return parser


def main(argv=None):
    """Write the STAR signal bedGraphs of a BAM as bigWigs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    print(args)
    star_return_code = call_star(args.bamfile, args.is_stranded, args.output_dir)

    try:
        assert star_return_code == 0
    except AssertionError:
        logger.exception("Building bedGraph had a problem, most likely out of memory.")
        sys.exit(1)

    kwargs = {
        "chrom_sizes": args.chrom_sizes,
        "threads": args.threads,
    }

    if args.is_stranded:
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.UniqueMultiple.str1.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_minusAll.bw"),
            **kwargs,
        )
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.Unique.str1.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_minusUniq.bw"),
            **kwargs,
        )
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.UniqueMultiple.str2.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_plusAll.bw"),
            **kwargs,
        )
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.Unique.str2.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_plusUniq.bw"),
            **kwargs,
        )
    else:
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.UniqueMultiple.str1.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_all.bw"),
            **kwargs,
        )
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.Unique.str1.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_uniq.bw"),
            **kwargs,
        )


def call_star(input_bam, is_stranded, output_dir):
    """Write the signal bedGraphs of a BAM with STAR.

    Parameters
    ----------
    input_bam : str
        Coordinate-sorted BAM.
    is_stranded : bool
        Write one bedGraph per strand.
    output_dir : str
        Directory for the bedGraphs.

    Returns
    -------
    int
        STAR exit code.
    """
    if is_stranded:
        strandedness = "stranded"
    else:
        strandedness = "unstranded"

    command = STAR_COMMAND.format(
        input_bam=input_bam,
        strandedness=strandedness.capitalize(),
        output_dir=output_dir,
    )
    logger.info("Running STAR command %s", command)
    return_code = subprocess.call(shlex.split(command))
    return return_code


def call_bg_to_bw(input_bg, out_fn, chrom_sizes, threads):
    """Sort a bedGraph, convert it to bigWig, and delete both bedGraphs.

    Parameters
    ----------
    input_bg : str
        The bedGraph from STAR.
    out_fn : str
        The bigWig to write.
    chrom_sizes : str
        Chromosome sizes file.
    threads : int
        Threads for `bigtools bedgraphtobigwig`.
    """
    sorted_bg = input_bg.replace(".bg", ".sorted.bg")

    bedgraph_cmd = f"bedtools sort -i {input_bg}"

    logger.info("Sorting bedgraph: %s", bedgraph_cmd)

    with open(sorted_bg, "w") as f:
        subprocess.call(shlex.split(bedgraph_cmd), stdout=f)
    command = f"bigtools bedgraphtobigwig --nthreads {threads} {sorted_bg} {chrom_sizes} {out_fn}"

    logger.info("Building bigWig: %s", command)

    subprocess.call(shlex.split(command))
    os.remove(input_bg)
    os.remove(sorted_bg)


if __name__ == "__main__":
    main()
