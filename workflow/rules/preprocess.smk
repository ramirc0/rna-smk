rule fastp:
    output:
        r1="results/processed/{sample}_R1.fastq.gz",
        r2="results/processed/{sample}_R2.fastq.gz",
        report="results/processed/{sample}.html",
        json="results/processed/{sample}.json"
    params:
        r1=lambda wildcards: READ_LOOKUP[wildcards.sample][0],
        r2=lambda wildcards: READ_LOOKUP[wildcards.sample][1],
    software: SOFTWARE_ENV
    log: "logs/processed/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Started fastp for sample: {wildcards.sample}"

        fastp -i {params.r1} -I {params.r2} \
              -o {output.r1} -O {output.r2} \
              -h {output.report} -j {output.json} \
              -w {resources.threads}
        
        echo "$(date): Completed fastp for sample: {wildcards.sample}"
        """
