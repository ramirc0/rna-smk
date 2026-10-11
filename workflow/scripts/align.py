"""Run the alignment (mapping) step of the ENCODE rna-seq-pipeline."""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import json
import logging
import os
import subprocess

from qc_utils import QCMetric
from qc_utils.parsers import parse_flagstats, parse_starlog

logger = logging.getLogger(__name__)

STAR_OPTIONS = [
    "--genomeLoad", "NoSharedMemory",
    "--outFilterMultimapNmax", "20",
    "--alignSJoverhangMin", "8",
    "--alignSJDBoverhangMin", "1",
    "--outFilterMismatchNmax", "999",
    "--outFilterMismatchNoverReadLmax", "0.04",
    "--alignIntronMin", "20",
    "--alignIntronMax", "1000000",
    "--alignMatesGapMax", "1000000",
    "--outSAMheaderHD", "@HD", "VN:1.4", "SO:coordinate",
    "--outSAMunmapped", "Within",
    "--outFilterType", "BySJout",
    "--outSAMattributes", "NH", "HI", "AS", "NM", "MD",
    "--outSAMtype", "BAM", "SortedByCoordinate",
    "--quantMode", "TranscriptomeSAM",
    "--sjdbScore", "1",
]  # fmt: skip


def run(command, **kwargs):
    """Log and run a command.

    Parameters
    ----------
    command : list of str
        Program and arguments.
    **kwargs
        Passed to `subprocess.run`.

    Returns
    -------
    subprocess.CompletedProcess
        The finished process.
    """
    logger.info("Running %s", " ".join(command))
    return subprocess.run(command, check=False, **kwargs)


def write_json(obj, path):
    """Write an object as JSON.

    Parameters
    ----------
    obj : object
        JSON-serializable object.
    path : str
        JSON file to write.
    """
    with open(path, "w") as f:
        json.dump(obj, f)


def build_parser():
    """Build the align.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the align.py flags.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", required=True, help="Output directory.")
    parser.add_argument("--fastq_r1", required=True, help="Gzipped read 1 FASTQ.")
    parser.add_argument("--fastq_r2", required=True, help="Gzipped read 2 FASTQ.")
    parser.add_argument("--indexdir", required=True, help="STAR index directory.")
    parser.add_argument(
        "--bamroot",
        required=True,
        help="Output name root. `foo` writes foo_genome.bam and foo_anno.bam.",
    )
    parser.add_argument("--ncpus", type=int, required=True, help="STAR threads.")
    parser.add_argument(
        "--ram_gb", type=int, required=True, help="RAM for BAM sorting, in GB."
    )
    return parser


def main(argv=None):
    """Align with STAR, then write the BAMs, flagstats and QC JSONs.

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
    prefix = os.path.join(out, args.bamroot)
    ncpus = str(args.ncpus)

    run(
        [
            "STAR",
            "--genomeDir", args.indexdir,
            "--readFilesIn", args.fastq_r1, args.fastq_r2,
            "--outFileNamePrefix", f"{out}/",
            "--readFilesCommand", "zcat",
            "--runThreadN", ncpus,
            *STAR_OPTIONS,
            "--limitBAMsortRAM", f"{args.ram_gb}000000000",
            "--outBAMsortingThreadN", ncpus,
            "--outBAMsortingBinsN", "100",
        ]
    )  # fmt: skip
    os.rename(f"{out}/Aligned.sortedByCoord.out.bam", f"{prefix}_genome.bam")
    os.rename(f"{out}/Log.final.out", f"{prefix}_Log.final.out")

    transcriptome_bam = f"{out}/Aligned.toTranscriptome.out.bam"
    validator = subprocess.run(
        ["rsem-sam-validator", transcriptome_bam],
        stdout=subprocess.PIPE,
        text=True,
        check=True,
    )
    # rsem-sam-validator exits 0 whether or not the BAM is valid, so the verdict
    # comes from its last line: "is valid!" or "is not valid!".
    if validator.stdout.strip().split("\n")[-1].endswith("is valid!"):
        logger.info("Transcriptome bam is already rsem-sorted.")
        os.rename(transcriptome_bam, f"{prefix}_anno.bam")
    else:
        logger.info("Transcriptome bam is not rsem-sorted.")
        run(["convert-sam-for-rsem", transcriptome_bam, f"{prefix}_anno"])

    for kind in ("genome", "anno"):
        flagstat = f"{prefix}_{kind}_flagstat.txt"
        with open(flagstat, "w") as f:
            run(["samtools", "flagstat", f"{prefix}_{kind}.bam"], stdout=f)
        qc = QCMetric(f"samtools_{kind}_flagstat", parse_flagstats(flagstat))
        write_json(qc.to_ordered_dict(), f"{prefix}_{kind}_flagstat.json")
    qc = QCMetric("star_log_qc", parse_starlog(f"{prefix}_Log.final.out"))
    write_json(qc.to_ordered_dict(), f"{prefix}_Log.final.json")


if __name__ == "__main__":
    main()
