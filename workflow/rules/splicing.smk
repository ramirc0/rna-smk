# leafcutter intron clusters: the inputs of differential splicing and sQTL mapping.
# Thresholds follow GTEx's sQTL pipeline (broadinstitute/gtex-pipeline, qtl/leafcutter).

LEAFCUTTER_DIR = f"{OUTDIR}/leafcutter"
INTRON_RATIOS = f"{LEAFCUTTER_DIR}/leafcutter_perind.counts.gz"

# STAR writes no XS tag here, so regtools takes each junction's strand from the
# library type. An unstranded library gives it none.
JUNCTION_STRAND = {"reverse": "RF", "forward": "FR"}
SPLICING_TARGETS = []
if config["sequencing"]["strand_direction"] in JUNCTION_STRAND:
    SPLICING_TARGETS += [
        f"{LEAFCUTTER_DIR}/leafcutter_perind_numers.counts.gz",
        EXON_TABLE,
    ]
    # One sample has no variation to normalize, and the phenotype script fails.
    if len(SAMPLES) > 1:
        SPLICING_TARGETS += [
            f"{INTRON_RATIOS}.qqnorm.bed.gz.tbi",
            f"{INTRON_RATIOS}.PCs",
        ]


rule junctions:
    input:
        f"{OUTDIR}/{{sample}}/align/{{sample}}_genome.bam",
    output:
        f"{OUTDIR}/{{sample}}/junctions/{{sample}}.junc.gz",
    params:
        strand=lambda w: JUNCTION_STRAND[config["sequencing"]["strand_direction"]],
    log:
        f"{LOGDIR}/junctions/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/junctions/{{sample}}.tsv",
    software:
        SPLICING_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT

        # regtools counts every alignment. MAPQ 255 is STAR's unique mapping.
        samtools view \
            -@ {threads:q} \
            -b \
            -q 255 \
            -o "$tmp/unique.bam" \
            {input:q}
        samtools index -@ {threads:q} "$tmp/unique.bam"

        regtools junctions extract \
            -a 8 \
            -m 50 \
            -M 500000 \
            -s {params.strand:q} \
            -o "$tmp/junctions.bed" \
            "$tmp/unique.bam"

        bgzip -c "$tmp/junctions.bed" > {output:q}

        sync {output:q}
        """


# A counts column is named after its junction file, up to ".junc".
rule intron_clusters:
    input:
        expand(f"{OUTDIR}/{{sample}}/junctions/{{sample}}.junc.gz", sample=SAMPLES),
    output:
        ratios=INTRON_RATIOS,
        counts=f"{LEAFCUTTER_DIR}/leafcutter_perind_numers.counts.gz",
        clusters=f"{LEAFCUTTER_DIR}/leafcutter_refined",
    log:
        f"{LOGDIR}/intron_clusters/all.txt",
    benchmark:
        f"{BENCHDIR}/intron_clusters/all.tsv"
    software:
        SPLICING_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # leafcutter writes a sorted table per sample, and merges in temp
        # directories it never removes.
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        export TMPDIR="$tmp"

        printf '%s\n' {input:q} > "$tmp/junction_files.txt"

        leafcutter-cluster \
            --juncfiles "$tmp/junction_files.txt" \
            --rundir "$tmp" \
            --outprefix leafcutter \
            --minclureads 30 \
            --mincluratio 0.001 \
            --maxintronlen 500000

        cp "$tmp/leafcutter_perind.counts.gz" {output.ratios:q}
        cp "$tmp/leafcutter_perind_numers.counts.gz" {output.counts:q}
        cp "$tmp/leafcutter_refined" {output.clusters:q}

        sync {output:q}
        """


rule splicing_phenotypes:
    input:
        INTRON_RATIOS,
    output:
        bed=f"{INTRON_RATIOS}.qqnorm.bed.gz",
        index=f"{INTRON_RATIOS}.qqnorm.bed.gz.tbi",
        pcs=f"{INTRON_RATIOS}.PCs",
    log:
        f"{LOGDIR}/splicing_phenotypes/all.txt",
    benchmark:
        f"{BENCHDIR}/splicing_phenotypes/all.tsv"
    software:
        SPLICING_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # The script writes next to its input.
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        ratios="$tmp/leafcutter_perind.counts.gz"
        ln -s "$(realpath {input:q})" "$ratios"

        # The script's own default never matches chr-prefixed names.
        echo chrY > "$tmp/excluded_chromosomes.txt"

        leafcutter-prepare-phenotype \
            --pcs 10 \
            --chromosome_blacklist "$tmp/excluded_chromosomes.txt" \
            "$ratios"

        # sort sizes its buffer from the node's RAM, not the job's.
        tables=("$ratios".qqnorm_*)
        {{
            head -n 1 "${{tables[0]}}"
            tail -q -n +2 "${{tables[@]}}" \
                | LC_ALL=C sort -t $'\t' -k1,1V -k2,2n -k3,3n -S 1G -T "$tmp"
        }} | bgzip > "$ratios.qqnorm.bed.gz"
        tabix -p bed "$ratios.qqnorm.bed.gz"

        cp "$ratios.qqnorm.bed.gz" {output.bed:q}
        cp "$ratios.qqnorm.bed.gz.tbi" {output.index:q}
        cp "$ratios.PCs" {output.pcs:q}

        sync {output:q}
        """


# Run-independent: labels leafcutter-ds clusters with gene names (--exon_file).
rule exon_table:
    input:
        REFERENCE_GTF,
    output:
        EXON_TABLE,
    log:
        f"{LOGDIR}/exon_table.txt",
    benchmark:
        f"{BENCHDIR}/exon_table.tsv"
    software:
        SPLICING_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        leafcutter-gtf-to-exons {input:q} {output:q}

        sync {output:q}
        """
