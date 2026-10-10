def get_bam_to_signals_rule():
    if config['sequencing']['is_stranded']:
        return {
            'output': {
                'minus_all': "results/signals/{sample}/{sample}_genome_minusAll.bw",
                'minus_unique': "results/signals/{sample}/{sample}_genome_minusUniq.bw",
                'plus_all': "results/signals/{sample}/{sample}_genome_plusAll.bw",
                'plus_unique': "results/signals/{sample}/{sample}_genome_plusUniq.bw",
            }
        }
    else:
        return {
            'output': {
                'unstranded_all': "results/signals/{sample}/{sample}_genome_all.bw",
                'unstranded_unique': "results/signals/{sample}/{sample}_genome_uniq.bw",
            }
        }

rule bam_to_signals:
    input:
        chromsizes=rules.index.output.chromsizes,
        genome_bam="results/align/{sample}/{sample}_genome.bam"
    output:
        **get_bam_to_signals_rule()['output'],
        signal="results/signals/{sample}/.continue",
    params:
        sample=lambda wildcards: wildcards.sample,
        is_stranded=config['sequencing']['is_stranded'],
        outdir=lambda wildcards: f"results/signals/{wildcards.sample}"
    software: SOFTWARE_ENV
    log: "logs/signals/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        ABS_LOG_PATH=$(realpath {log})
        echo "$(date): Running bam to singal conversion"
        
        samtools sort -@ {resources.threads} -m 1G -o {params.outdir}/temp_sorted.bam {input.genome_bam}

        python workflow/scripts/bam_to_signals.py \
            --threads {resources.threads} \
            --output_dir {params.outdir} \
            --bamfile {params.outdir}/temp_sorted.bam \
            --chrom_sizes {input.chromsizes} \
            --is_stranded {params.is_stranded} \
            --bamroot {params.sample}_genome
        
        rm {params.outdir}/temp_sorted.bam

        echo "$(date): Finished bam to singal conversion"
        touch {output.signal}
        """