# Genome BAM to bigWig signal tracks.


rule bam_to_signals:
    input:
        chromsizes=f"{STAR_INDEX}/chrNameLength.txt",
        genome_bam=f"{OUTDIR}/{{sample}}/align/{{sample}}_genome.bam",
    output:
        SIGNAL_TRACKS,
    params:
        outdir=lambda w: f"{OUTDIR}/{w.sample}/signals",
        strandedness=lambda w: (
            "stranded" if config["sequencing"]["is_stranded"] else "unstranded"
        ),
    log:
        f"{LOGDIR}/bam_to_signals/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/bam_to_signals/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/bam_to_signals.py \
            --threads {threads:q} \
            --output_dir {params.outdir:q} \
            --bamfile {input.genome_bam:q} \
            --chrom_sizes {input.chromsizes:q} \
            --strandedness {params.strandedness:q} \
            --bamroot {wildcards.sample:q}_genome
        """
