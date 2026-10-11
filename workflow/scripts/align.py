"""Run the alignment (mapping) step of the ENCODE rna-seq-pipeline."""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import json
import logging
import os
import pathlib
import random
import re
import shlex
import shutil
import subprocess
from abc import ABC, abstractmethod

from qc_utils import QCMetric
from qc_utils.parsers import parse_flagstats, parse_starlog

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)


def choices(population, k):
    """Draw `k` elements from `population` with replacement.

    Parameters
    ----------
    population : sequence
        Elements to draw from.
    k : int
        Number of draws.

    Returns
    -------
    list
        The drawn elements.
    """
    result = []
    for i in range(k):
        result.append(random.choice(population))
    return result


def get_tmp_file_name(extension=".fastq.gz"):
    """Return an unused file name of 20 random lowercase letters.

    Parameters
    ----------
    extension : str, default ".fastq.gz"
        Suffix of the name.

    Returns
    -------
    str
        A name that no file in the working directory has.
    """
    while True:
        tmp_name = "".join(choices("abcdefghijklmnopqrstuvwxyz", k=20)) + extension
        if pathlib.Path(tmp_name).exists():
            continue
        else:
            return tmp_name


def concatenate_files(input_files):
    """Concatenate files into a new file in the working directory.

    Parameters
    ----------
    input_files : list of str
        Files to concatenate, in order.

    Returns
    -------
    str
        Name of the concatenated file.
    """
    result_filename = get_tmp_file_name()
    with open(result_filename, "wb") as out_fp:
        logger.info(f"merging files into {result_filename!r}")
        for fastq in input_files:
            with open(fastq, "rb") as add_on:
                logger.info(f"merging {fastq!r} next")
                shutil.copyfileobj(add_on, out_fp)
                logger.info(f"merging {fastq!r} success")
    logger.info(f"merge complete, result is in {result_filename!r}")
    return result_filename


def get_flagstats(input_path, output_path):
    """Write the `samtools flagstat` report of a BAM.

    Parameters
    ----------
    input_path : str
        BAM to summarize.
    output_path : str
        Text report to write.
    """
    command = f"samtools flagstat {input_path}"
    logger.info("Getting samtools flagstats for %s", input_path)
    process = subprocess.run(shlex.split(command), stdout=subprocess.PIPE, check=False)
    with open(output_path, "w") as f:
        f.write(process.stdout.decode())


def write_json(input_obj, output_path):
    """Write an object as JSON.

    Parameters
    ----------
    input_obj : object
        JSON-serializable object.
    output_path : str
        JSON file to write.
    """
    with open(output_path, "w") as fp:
        json.dump(input_obj, fp)


class StarAligner(ABC):
    """Settings and run step of a STAR job.

    Parameters
    ----------
    ncpus : int
        STAR threads.
    ram_gb : int
        RAM for BAM sorting, in GB.
    indexdir : str
        STAR index directory.
    outdir : str
        Output directory.
    """

    def __init__(self, ncpus, ram_gb, indexdir, outdir):
        """Store the STAR settings.

        Parameters
        ----------
        ncpus : int
            STAR threads.
        ram_gb : int
            RAM for BAM sorting, in GB.
        indexdir : str
            STAR index directory.
        outdir : str
            Output directory.
        """
        self.ncpus = ncpus
        self.ram_gb = ram_gb
        self.indexdir = indexdir
        self.outdir = outdir

    def run(self):
        """Run the STAR command."""
        logger.info("running STAR command %s", " ".join(self.command))
        subprocess.call(self.command)

    @property
    @abstractmethod
    def command_string(self):
        """Return the STAR command template.

        Returns
        -------
        str
            Command with `str.format` fields for the settings.
        """

    @abstractmethod
    def format_command_string(self):
        """Fill the STAR command template.

        Returns
        -------
        str
            The STAR command.
        """


class PairedEndStarAligner(StarAligner):
    """STAR job for a paired-end library.

    Parameters
    ----------
    fastqs : list of str
        Read 1 and read 2 FASTQs.
    ncpus : int
        STAR threads.
    ram_gb : int
        RAM for BAM sorting, in GB.
    indexdir : str
        STAR index directory.
    outdir : str
        Output directory.
    """

    command_string = """STAR --genomeDir {indexdir} \
    --readFilesIn {read1_fq_gz} {read2_fq_gz} \
    --outFileNamePrefix {outdir}/ \
    --readFilesCommand zcat \
    --runThreadN {ncpus} \
    --genomeLoad NoSharedMemory \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMheaderCommentFile COfile.txt \
    --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
    --sjdbScore 1 \
    --limitBAMsortRAM {ram_gb}000000000 \
    --outBAMsortingThreadN {ncpus} \
    --outBAMsortingBinsN 100"""

    def __init__(self, fastqs, ncpus, ram_gb, indexdir, outdir):
        """Build the STAR command for a FASTQ pair.

        Parameters
        ----------
        fastqs : list of str
            Read 1 and read 2 FASTQs.
        ncpus : int
            STAR threads.
        ram_gb : int
            RAM for BAM sorting, in GB.
        indexdir : str
            STAR index directory.
        outdir : str
            Output directory.
        """
        super().__init__(ncpus, ram_gb, indexdir, outdir)
        self.fastq_read1 = fastqs[0]
        self.fastq_read2 = fastqs[1]
        self.command = shlex.split(
            self.format_command_string(type(self).command_string)
        )

    def format_command_string(self, input_string):
        """Fill a STAR command template with this job's settings.

        Parameters
        ----------
        input_string : str
            Command template.

        Returns
        -------
        str
            The STAR command.
        """
        logger.info(
            f"Formatting with arguments: \nR1: {self.fastq_read1} \nR2: {self.fastq_read2} \nNCPUS: {self.ncpus} \nRAMGB: {self.ram_gb} \nINDEXDIR: {self.indexdir} \nOUTDIR: {self.outdir}"
        )

        cmd = input_string.format(
            read1_fq_gz=self.fastq_read1,
            read2_fq_gz=self.fastq_read2,
            ncpus=self.ncpus,
            ram_gb=self.ram_gb,
            indexdir=self.indexdir,
            outdir=self.outdir,
        )
        return cmd


def build_parser():
    """Build the align.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the align.py flags.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", help="Path to output directory.", required=True)
    parser.add_argument(
        "--fastqs_r1", nargs="+", help="Input gzipped fastq(s) belonging to read1"
    )
    parser.add_argument(
        "--fastqs_r2", nargs="+", help="Input gzipped fastq(s) belonging to read2"
    )
    parser.add_argument(
        "--indexdir",
        type=str,
        help="Directory of the STAR index (genomeDir).",
        default="out",
    )
    parser.add_argument(
        "--bamroot",
        type=str,
        help="""
             Root name for output bams. For example out_bam
             will create out_bam_genome.bam and out_bam_anno.bam
             """,
        default="out_bam",
    )
    parser.add_argument(
        "--ncpus", type=int, help="Number of cpus available.", default=4
    )
    parser.add_argument(
        "--ram_gb", type=int, help="Amount of RAM available in GB.", default=8
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
    fastqs = [
        concatenate_files(reads) if len(reads) > 1 else reads[0]
        for reads in (args.fastqs_r1, args.fastqs_r2)
    ]
    aligner = PairedEndStarAligner(
        fastqs, args.ncpus, args.ram_gb, args.indexdir, args.output_dir
    )
    aligner.run()
    output_dir = args.output_dir
    genome_bam_path = os.path.join(output_dir, args.bamroot + "_genome.bam")
    anno_bam_path = os.path.join(output_dir, args.bamroot + "_anno.bam")
    genome_flagstat_path = os.path.join(
        output_dir, args.bamroot + "_genome_flagstat.txt"
    )
    anno_flagstat_path = os.path.join(output_dir, args.bamroot + "_anno_flagstat.txt")
    star_log_path = os.path.join(output_dir, args.bamroot + "_Log.final.out")
    os.rename(
        os.path.join(output_dir, "Aligned.sortedByCoord.out.bam"), genome_bam_path
    )
    os.rename(os.path.join(output_dir, "Log.final.out"), star_log_path)
    transcriptome_bam = os.path.join(output_dir, "Aligned.toTranscriptome.out.bam")
    rsem_check_cmd = f"rsem-sam-validator {transcriptome_bam}"
    rsem_output = subprocess.check_output(shlex.split(rsem_check_cmd))
    # rsem-sam-validator exits 0 whether or not the BAM is valid, so the verdict
    # comes from its last line: "is valid!" or "is not valid!".
    rsem_valid = rsem_output.decode().strip().split("\n")[-1].endswith("is valid!")
    if rsem_valid:
        logger.info("Transcriptome bam is already rsem-sorted.")
        os.rename(
            os.path.join(output_dir, "Aligned.toTranscriptome.out.bam"), anno_bam_path
        )
    else:
        logger.info("Transcriptome bam is not rsem-sorted.")
        rsem_sort_cmd = "convert-sam-for-rsem {input} {output}".format(
            input=os.path.join(output_dir, "Aligned.toTranscriptome.out.bam"),
            output=anno_bam_path.rstrip(".bam"),
        )
        logger.info("Running %s", rsem_sort_cmd)
        subprocess.call(shlex.split(rsem_sort_cmd))
    get_flagstats(genome_bam_path, genome_flagstat_path)
    get_flagstats(anno_bam_path, anno_flagstat_path)
    anno_flagstat_content = parse_flagstats(anno_flagstat_path)
    genome_flagstat_content = parse_flagstats(genome_flagstat_path)
    star_log_content = parse_starlog(star_log_path)
    anno_flagstat_qc = QCMetric("samtools_anno_flagstat", anno_flagstat_content)
    genome_flagstat_qc = QCMetric("samtools_genome_flagstat", genome_flagstat_content)
    star_log_qc = QCMetric("star_log_qc", star_log_content)
    write_json(
        anno_flagstat_qc.to_ordered_dict(),
        re.sub(r"\.txt$", ".json", anno_flagstat_path),
    )
    write_json(
        genome_flagstat_qc.to_ordered_dict(),
        re.sub(r"\.txt$", ".json", genome_flagstat_path),
    )
    write_json(star_log_qc.to_ordered_dict(), re.sub(r"\.out$", ".json", star_log_path))


if __name__ == "__main__":
    main()
