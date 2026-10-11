"""Run the alignment (mapping) step of the ENCODE rna-seq-pipeline."""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import json
import logging
import os
import re
import subprocess

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
    """Log and run a command, raising on a nonzero exit.

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
    return subprocess.run(command, check=True, **kwargs)


# JSON key of each `samtools flagstat` line that ENCODE's QC reports. ENCODE's
# parser let "primary duplicates" overwrite "duplicates", so that line wins.
FLAGSTAT_KEYS = {
    "in total (QC-passed reads + QC-failed reads)": "total",
    "primary duplicates": "duplicates",
    "mapped": "mapped",
    "paired in sequencing": "paired",
    "read1": "read1",
    "read2": "read2",
    "properly paired": "paired_properly",
    "with itself and mate mapped": "with_itself",
    "singletons": "singletons",
    "with mate mapped to a different chr": "diff_chroms",
}


def parse_flagstat(path):
    """Read the counts and percentages of a `samtools flagstat` report.

    Parameters
    ----------
    path : str
        Text report of `samtools flagstat`.

    Returns
    -------
    dict of str to int or float
        Counts, QC-failed counts and percentages, keyed as in `FLAGSTAT_KEYS`.
    """
    metrics = {}
    with open(path) as f:
        for line in f:
            m = re.fullmatch(r"(\d+) \+ (\d+) (.+?)(?: \(([\d.]+)%.*\))?", line.strip())
            key = m and FLAGSTAT_KEYS.get(m[3])
            if key:
                metrics[key] = int(m[1])
                metrics[f"{key}_qc_failed"] = int(m[2])
                if m[4]:
                    metrics[f"{key}_pct"] = float(m[4])
    return metrics


def parse_star_log(path):
    """Read the `name | value` lines of a STAR `Log.final.out`.

    Parameters
    ----------
    path : str
        STAR `Log.final.out`.

    Returns
    -------
    dict of str to str
        Each metric name and its value, both stripped.
    """
    with open(path) as f:
        rows = [line.split("|") for line in f if "|" in line]
    return {row[0].strip(): row[1].strip() for row in rows}


def write_qc_json(name, metrics, path):
    """Write metrics as ENCODE's `{name: {metric: value}}` QC JSON.

    Parameters
    ----------
    name : str
        Top-level key.
    metrics : dict
        Metrics, written with sorted keys.
    path : str
        JSON file to write.
    """
    with open(path, "w") as f:
        json.dump({name: metrics}, f, sort_keys=True)


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
            run(["samtools", "flagstat", "-@", ncpus, f"{prefix}_{kind}.bam"], stdout=f)
        write_qc_json(
            f"samtools_{kind}_flagstat",
            parse_flagstat(flagstat),
            f"{prefix}_{kind}_flagstat.json",
        )
    write_qc_json(
        "star_log_qc",
        parse_star_log(f"{prefix}_Log.final.out"),
        f"{prefix}_Log.final.json",
    )


if __name__ == "__main__":
    main()
