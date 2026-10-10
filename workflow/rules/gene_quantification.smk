# RSEM reference and gene and isoform quantification.


# Run-independent: built once under results/refs/, reused across runs.
rule rsem_index:
    input:
        fasta=REFERENCE_FASTA,
        gtf=REFERENCE_GTF,
    output:
        RSEM_INDEX_FILES,
    params:
        prefix=RSEM_INDEX,
    log:
        f"{LOGDIR}/rsem_index.txt",
    benchmark:
        f"{BENCHDIR}/rsem_index.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        ulimit -n 65536

        rsem-prepare-reference \
            --gtf {input.gtf:q} \
            -p {threads:q} \
            {input.fasta:q} \
            {params.prefix:q}
        """


# RSEM's Gibbs sampling and credibility intervals depend on -p.
rule rsem_quant:
    input:
        index=RSEM_INDEX_FILES,
        anno_bam=f"{OUTDIR}/{{sample}}/align/{{sample}}_anno.bam",
    output:
        genes=f"{OUTDIR}/{{sample}}/rsem_quant/{{sample}}.genes.results",
        isoforms=f"{OUTDIR}/{{sample}}/rsem_quant/{{sample}}.isoforms.results",
        genes_detected=f"{OUTDIR}/{{sample}}/rsem_quant/{{sample}}_no_genes_detected.json",
    params:
        index=RSEM_INDEX,
        outprefix=lambda w: f"{OUTDIR}/{w.sample}/rsem_quant/{w.sample}",
        endedness=lambda w: config["sequencing"]["endedness"],
        strand_direction=lambda w: config["sequencing"]["strand_direction"],
        rsem_seed=lambda w: config["sequencing"]["rsem_seed"],
    log:
        f"{LOGDIR}/rsem_quant/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/rsem_quant/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        ulimit -n 65536

        python workflow/scripts/rsem.py \
            --rsem_index {params.index:q} \
            --anno_bam {input.anno_bam:q} \
            --endedness {params.endedness:q} \
            --read_strand {params.strand_direction:q} \
            --rnd_seed {params.rsem_seed:q} \
            --ncpus {threads:q} \
            --mem_mb {resources.mem_mb:q} \
            --outprefix {params.outprefix:q}
        """
