from pathlib import Path

rule rsem_index:
    output:
        signal="results/rsem_index.continue",
        index_dir=directory("results/rsem_index")
    params:
        fasta=config['input_files']['reference_fasta'],
        gtf=config['input_files']['gtf_file'],
        transcript_to_gene_map=(f"--transcript-to-gene-map {Path(config['input_files']['knownIsoforms']).name}", config['input_files']['knownIsoforms']) if config['input_files'].get('knownIsoforms') else "",
        rsem_index_prefix="results/rsem_index/rsem"
    software: SOFTWARE_ENV
    log: "logs/rsem_index.log"
    shell:
        """
        exec >> {log} 2>&1

        echo "$(date): Creating rsem index"
        ulimit -n 65536
        echo "$(date): ulimit set to $(ulimit -n)"

        mkdir -p results/rsem_index

        rsem-prepare-reference \
            --gtf {params.gtf} \
            -p {resources.threads} \
            {params.fasta} \
            {params.rsem_index_prefix}

        echo "$(date): Finished creating rsem index"
        touch {output.signal}
        """

rule rsem_quant:
    input:
        signal=rules.rsem_index.output.signal,
        index_dir=rules.rsem_index.output.index_dir,
        anno_bam="results/align/{sample}/{sample}_anno.bam",
    output:
        signal="results/rsem_quant/{sample}.continue",
    params:
        outprefix=lambda wildcards: f"results/rsem_quant/{wildcards.sample}",
        endedness=config['sequencing']['endedness'],
        strand_direction=config['sequencing']['strand_direction'],
        rsem_seed=config['sequencing']['rsem_seed']
    software: SOFTWARE_ENV
    log: "logs/rsem_quant/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        
        echo "$(date): Running rsem quantification"
        ulimit -n 65536
        echo "$(date): ulimit set to $(ulimit -n)"

        python workflow/scripts/rsem.py \
            --rsem_index {input.index_dir}/rsem \
            --anno_bam {input.anno_bam} \
            --endedness {params.endedness} \
            --read_strand {params.strand_direction} \
            --rnd_seed {params.rsem_seed} \
            --ncpus {resources.threads} \
            --mem_mb {resources.mem_mb} \
            --outprefix {params.outprefix}

        echo "$(date): Finished rsem quantification"
        touch {output.signal}
        """
