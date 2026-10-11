# Strandedness inference.


# Run-independent: built once under results/refs/, reused across runs.
rule gene_model:
    input:
        REFERENCE_GTF,
    output:
        GENE_MODEL,
    log:
        f"{LOGDIR}/gene_model.txt",
    benchmark:
        f"{BENCHDIR}/gene_model.tsv",
    software:
        QC_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gtfToGenePred -genePredExt {input:q} stdout \
            | genePredToBed stdin {output:q}
        """


rule strandedness:
    input:
        bam=f"{OUTDIR}/{{sample}}/align/{{sample}}_genome.bam",
        gene_model=GENE_MODEL,
    output:
        f"{OUTDIR}/{{sample}}/strandedness/{{sample}}.txt",
    log:
        f"{LOGDIR}/strandedness/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/strandedness/{{sample}}.tsv",
    software:
        QC_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/get_strandedness.py \
            --bam {input.bam:q} \
            --gene_model {input.gene_model:q} \
            --output {output:q}
        """

