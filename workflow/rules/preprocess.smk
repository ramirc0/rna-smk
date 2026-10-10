# Adapter and quality trimming.


rule fastp:
    input:
        r1=lambda w: R1_OF[w.sample],
        r2=lambda w: R2_OF[w.sample],
    output:
        r1=f"{OUTDIR}/{{sample}}/processed/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{{sample}}/processed/{{sample}}_R2.fastq.gz",
        html=f"{OUTDIR}/{{sample}}/processed/{{sample}}.html",
        json=f"{OUTDIR}/{{sample}}/processed/{{sample}}.json",
    log:
        f"{LOGDIR}/fastp/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/fastp/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        fastp \
            -i {input.r1:q} \
            -I {input.r2:q} \
            -o {output.r1:q} \
            -O {output.r2:q} \
            -h {output.html:q} \
            -j {output.json:q} \
            -w {threads:q}
        """
