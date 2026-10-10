# Strandedness inference and MAD QC between replicates.


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


rule mad_qc:
    input:
        tsv1=f"{OUTDIR}/{{sample}}/rsem_quant/{{sample}}.genes.results",
        tsv2=f"{OUTDIR}/{{other}}/rsem_quant/{{other}}.genes.results",
    output:
        json=f"{OUTDIR}/{{sample}}/mad_qc/{{other}}.json",
        plot=f"{OUTDIR}/{{sample}}/mad_qc/{{other}}.svg",
        plot_png=f"{OUTDIR}/{{sample}}/mad_qc/{{other}}.png",
    log:
        f"{LOGDIR}/mad_qc/{{sample}}/{{other}}.txt",
    benchmark:
        f"{BENCHDIR}/mad_qc/{{sample}}/{{other}}.tsv",
    software:
        QC_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/mad_qc.py \
            --tsv1 {input.tsv1:q} \
            --tsv2 {input.tsv2:q} \
            --json {output.json:q} \
            --plot {output.plot:q}
        """
