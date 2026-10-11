# STAR genome index and alignment.


# Run-independent: built once under results/refs/, reused across runs.
rule index:
    input:
        fasta=REFERENCE_FASTA,
        gtf=REFERENCE_GTF,
    output:
        index_dir=directory(STAR_INDEX),
        chromsizes=f"{STAR_INDEX}/chrNameLength.txt",
    params:
        sjdb_overhang=lambda w: config["sequencing"]["sjdb_overhang"],
    log:
        f"{LOGDIR}/index.txt",
    benchmark:
        f"{BENCHDIR}/index.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        ulimit -n 65536

        STAR \
            --runThreadN {threads:q} \
            --runMode genomeGenerate \
            --genomeDir {output.index_dir:q} \
            --genomeFastaFiles {input.fasta:q} \
            --sjdbGTFfile {input.gtf:q} \
            --sjdbOverhang {params.sjdb_overhang:q}
        """


rule align:
    input:
        index_dir=STAR_INDEX,
        r1=f"{OUTDIR}/{{sample}}/processed/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{{sample}}/processed/{{sample}}_R2.fastq.gz",
    output:
        genome_bam=f"{OUTDIR}/{{sample}}/align/{{sample}}_genome.bam",
        anno_bam=f"{OUTDIR}/{{sample}}/align/{{sample}}_anno.bam",
        genome_flagstat=f"{OUTDIR}/{{sample}}/align/{{sample}}_genome_flagstat.txt",
        anno_flagstat=f"{OUTDIR}/{{sample}}/align/{{sample}}_anno_flagstat.txt",
        star_log=f"{OUTDIR}/{{sample}}/align/{{sample}}_Log.final.out",
        genome_flagstat_json=f"{OUTDIR}/{{sample}}/align/{{sample}}_genome_flagstat.json",
        anno_flagstat_json=f"{OUTDIR}/{{sample}}/align/{{sample}}_anno_flagstat.json",
        star_log_json=f"{OUTDIR}/{{sample}}/align/{{sample}}_Log.final.json",
    params:
        outdir=lambda w: f"{OUTDIR}/{w.sample}/align",
        ram_gb=lambda w, resources: resources.mem_mb // 1000,
    log:
        f"{LOGDIR}/align/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/align/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        ulimit -n 65536

        python workflow/scripts/align.py \
            --output_dir {params.outdir:q} \
            --fastq_r1 {input.r1:q} \
            --fastq_r2 {input.r2:q} \
            --bamroot {wildcards.sample:q} \
            --indexdir {input.index_dir:q} \
            --ram_gb {params.ram_gb:q} \
            --ncpus {threads:q}

        samtools quickcheck {output.genome_bam:q} {output.anno_bam:q} && \
            echo "$(date): samtools quickcheck OK" || echo "$(date): samtools quickcheck FAIL"
        """
