ATTEMPT_TO_QUEUE = {
    1: ("12hours", 720),
    2: ("5days", 7200),
    3: ("long", 43200)
}

def get_slurm_partition(wildcards, attempt):
    return ATTEMPT_TO_QUEUE[attempt][0]

def get_runtime(wildcards, attempt):
    return ATTEMPT_TO_QUEUE[attempt][1]

rule index:
    output:
        chromsizes="results/index/chrNameLength.txt",
        signal="results/index.continue",
        index_dir=directory("results/index/"),
    params:
        sample=FIRST_SAMPLE,
        fasta=config['input_files']['reference_fasta'],
        gtf=config['input_files']['gtf_file'],
    software: SOFTWARE_ENV
    log: "logs/index.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Creating STAR index"
        ulimit -n 65536
        echo "$(date): ulimit set to $(ulimit -n)"

        AVG_READ_LENGTH=$(seqkit stats {params.sample} | awk 'NR>1 {{print int($6)}}')

        STAR \
            --runThreadN {resources.threads} \
            --runMode genomeGenerate \
            --genomeDir {output.index_dir} \
            --genomeFastaFiles {params.fasta} \
            --sjdbGTFfile {params.gtf} \
            --sjdbOverhang $(($AVG_READ_LENGTH - 1))

        echo "$(date): Finished creating STAR index"
        touch {output.signal}
        """

rule align:
    input:
        signal=rules.index.output.signal,
        index_dir=rules.index.output.index_dir,
        r1="results/processed/{sample}_R1.fastq.gz",
        r2="results/processed/{sample}_R2.fastq.gz"
    output:
        genome_bam="results/align/{sample}/{sample}_genome.bam",
        anno_bam="results/align/{sample}/{sample}_anno.bam",
        genome_flagstat="results/align/{sample}/{sample}_genome_flagstat.txt",
        anno_flagstat="results/align/{sample}/{sample}_anno_flagstat.txt",
        log="results/align/{sample}/{sample}_Log.final.out",
        genome_flagstat_json="results/align/{sample}/{sample}_genome_flagstat.json",
        anno_flagstat_json="results/align/{sample}/{sample}_anno_flagstat.json",
        log_json="results/align/{sample}/{sample}_Log.final.json",
        outdir=directory("results/align/{sample}/"),
        signal="results/align/{sample}/.continue"
    retries: 3
    params:
        endedness=config['sequencing']['endedness'],
    resources:
        slurm_partition=get_slurm_partition,
        runtime=get_runtime,
    software: SOFTWARE_ENV
    log: "logs/align/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Running STAR alignment"
        ulimit -n 65536
        echo "$(date): ulimit set to $(ulimit -n)"

        python workflow/scripts/align.py \
            --output_dir {output.outdir} \
            --fastqs_r1 {input.r1} \
            --fastqs_r2 {input.r2} \
            --bamroot {wildcards.sample} \
            --indexdir {input.index_dir} \
            --endedness {params.endedness} \
            --ram_gb {resources.mem_gb} \
            --ncpus {resources.threads} 
        
        samtools quickcheck {output.genome_bam} {output.anno_bam} && \
            echo "$(date): samtools quickcheck OK" || echo "$(date): samtools quickcheck FAIL"

        echo "$(date): Finished STAR alignment"
        touch {output.signal}
        """