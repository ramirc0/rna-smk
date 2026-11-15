# The run's QC report: every sample's QC in one table, graded on ENCODE's bulk
# RNA-seq portal audits (workflow/report/qc.yaml), with interactive figures.

import subprocess

REPORT_DIR = f"{OUTDIR}/report"
# The header names the pipeline commit, so a new commit renders the page again.
PIPELINE_VERSION = subprocess.run(
    ["git", "-C", workflow.basedir, "describe", "--always", "--dirty"],
    capture_output=True,
    text=True,
).stdout.strip()


# Run-independent: the gene types and chromosomes behind the composition metrics.
rule gene_table:
    input:
        REFERENCE_GTF,
    output:
        GENE_TABLE,
    log:
        f"{LOGDIR}/gene_table.txt",
    benchmark:
        f"{BENCHDIR}/gene_table.tsv"
    software:
        REPORT_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/gene_table.py \
            --gtf {input:q} \
            --output {output:q}

        sync {output:q}
        """


rule gather_qc:
    input:
        gene_table=GENE_TABLE,
        qc=expand(
            [
                f"{OUTDIR}/{{sample}}/processed/{{sample}}.json",
                f"{OUTDIR}/{{sample}}/align/{{sample}}_Log.final.json",
                f"{OUTDIR}/{{sample}}/align/{{sample}}_genome_flagstat.json",
                f"{OUTDIR}/{{sample}}/rsem_quant/{{sample}}.genes.results",
                f"{OUTDIR}/{{sample}}/rsem_quant/{{sample}}_no_genes_detected.json",
                f"{OUTDIR}/{{sample}}/strandedness/{{sample}}_fractions.tsv",
            ],
            sample=SAMPLES,
        ),
    output:
        f"{REPORT_DIR}/qc_metrics.tsv",
    params:
        run=OUTDIR,
        samples=SAMPLES,
    log:
        f"{LOGDIR}/gather_qc/all.txt",
    benchmark:
        f"{BENCHDIR}/gather_qc/all.tsv"
    software:
        REPORT_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/gather_qc.py \
            --run-dir {params.run:q} \
            --samples {params.samples:q} \
            --gene-table {input.gene_table:q} \
            --threads {threads:q} \
            --output {output:q}

        sync {output:q}
        """


rule qc_figures:
    input:
        gene_table=GENE_TABLE,
        qc=expand(
            [
                f"{OUTDIR}/{{sample}}/processed/{{sample}}.json",
                f"{OUTDIR}/{{sample}}/align/{{sample}}_Log.final.json",
                f"{OUTDIR}/{{sample}}/rsem_quant/{{sample}}.genes.results",
                f"{OUTDIR}/{{sample}}/strandedness/{{sample}}_fractions.tsv",
            ],
            sample=SAMPLES,
        ),
    output:
        expand(f"{REPORT_DIR}/samples/{{sample}}.js", sample=SAMPLES),
    params:
        run=OUTDIR,
        report_dir=REPORT_DIR,
        samples=SAMPLES,
    log:
        f"{LOGDIR}/qc_figures/all.txt",
    benchmark:
        f"{BENCHDIR}/qc_figures/all.tsv"
    software:
        REPORT_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/qc_figures.py \
            --run-dir {params.run:q} \
            --report-dir {params.report_dir:q} \
            --samples {params.samples:q} \
            --gene-table {input.gene_table:q} \
            --threads {threads:q}

        sync {output:q}
        """


rule qc_report:
    input:
        spec=f"{workflow.basedir}/report/qc.yaml",
        metrics=f"{REPORT_DIR}/qc_metrics.tsv",
        figures=expand(f"{REPORT_DIR}/samples/{{sample}}.js", sample=SAMPLES),
    output:
        f"{REPORT_DIR}/qc_report.html",
    params:
        run=OUTDIR,
        run_id=RUN_ID,
        version=f"rna-smk {PIPELINE_VERSION}",
    log:
        f"{LOGDIR}/qc_report/all.txt",
    benchmark:
        f"{BENCHDIR}/qc_report/all.tsv"
    software:
        REPORT_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        qc-report \
            --spec {input.spec:q} \
            --metrics {input.metrics:q} \
            --run-dir {params.run:q} \
            --run-id {params.run_id:q} \
            --version {params.version:q} \
            --output {output:q}

        sync {output:q}
        """
