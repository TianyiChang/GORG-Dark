# produce: the numbers of mapped and total (rarefied not original), and post-QC reads for each sample

import pandas as pd

metag_full_list = pd.read_table('../../metadata/local_metag_list.txt')
metag_ids = list(metag_full_list.ids)

# sra_run_ids, = glob_wildcards("../../metag_wildcards/from_sra/{sra}.wildcard")
# metag_ids, = glob_wildcards("../../metag_wildcards/{metag}.wildcard")

sag_id, = glob_wildcards("../../low_entropy/test/gd_sags/{sag}.fasta")
bac_id, = glob_wildcards("../../low_entropy/sags/bac/{bac}.fasta")
arc_id, = glob_wildcards("../../low_entropy/sags/arc/{arc}.fasta")

ref_genomes = [
    'outside_gorg_tropics', 'outside_omd_mags', 'outside_acinas_2020',
    'gorg_v4_concat', 'gorg_v4_omd', 'omd_m_dark', 'gorg_v4_omd_dark']

rule all:
    input:
        expand("../GRCh38_mapped_seq/{metag}.fa", metag=metag_ids),
        expand("../{ref_genomes}/stats/n_total_read.csv", ref_genomes=ref_genomes),
        expand("../{ref_genomes}/stats/n_mapped_read.csv", ref_genomes=ref_genomes),
        expand("../{ref_genomes}/post_qc_summary/metag_post_qc_read.csv",
             ref_genomes=ref_genomes)

rule trimmomatic:
    input:
        "../../metag/from_collab/{metag}_1.fastq.gz",
        "../../metag/from_collab/{metag}_2.fastq.gz"
    output:
        temp("../trimmomatic/{metag}_1_paired.fq.gz"),
        temp("../trimmomatic/{metag}_1_unpaired.fq.gz"),
        temp("../trimmomatic/{metag}_2_paired.fq.gz"),
        temp("../trimmomatic/{metag}_2_unpaired.fq.gz")
    conda:
        "trimmomatic"
    threads: 12
    params: mem="5g"
    shell:
        """
        trimmomatic PE -phred33 -threads 12 \
        {input[0]} {input[1]} \
        {output[0]} {output[1]} {output[2]} {output[3]} \
        ILLUMINACLIP:/home/tchang/.conda/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

rule read_merge:
    input:
        f="../trimmomatic/{metag}_1_paired.fq.gz",
        r="../trimmomatic/{metag}_2_paired.fq.gz"
    output:
        temp("../merged/{metag}.extendedFrags.fastq.gzip")
    conda: "../envs/flash.yml"
    threads: 1
    params:
        outdir="../merged",
        mem="5g"
    shell:
        """
        (flash -x 0.05 -m 20 -M 999 \
        {input.f} {input.r} -d {params.outdir} \
        -o {wildcards.metag} --compress-prog=gzip &&
        rm {params.outdir}/{wildcards.metag}.notCombined*) ||
            touch {params.outdir}/{wildcards.metag}_all_reads_joined
        """

rule read_rarefy:
    input:
        "../merged/{metag}.extendedFrags.fastq.gzip"
    output:
        "../rarefy/{metag}.fq"
    params:
        seed="940216",
        n_read="1000000",
        outdir="../rarefy",
        mem="5g"
    conda: "../envs/seqtk.yml"
    threads: 1
    shell:
        """
        if [[ $(gunzip -c {input} | wc -l) -ge 4000000 ]]; then
            seqtk sample -s{params.seed} {input} \
            {params.n_read} > {output}
        else
            cp {input} {params.outdir}/{wildcards.metag}.fq.gz
            gunzip {params.outdir}/{wildcards.metag}.fq.gz
        fi
        """

rule rm_lc_read:
    input:
        rscript="../../../scripts/pcoa_sag_abund/frag_recruit_lcr_reads.r",
        fq="../rarefy/{metag}.fq"
    output:
        temp("../lc_rmd_reads/{metag}_lc_rmd.fq")
    threads: 4
    params:
        distinct_2mer="0.15",
        base_freq="0.02",
        mem="1g"
    shell:
        """
        {input.rscript} \
            {threads} {input.fq} {params.distinct_2mer} \
            {params.base_freq} {output}
        """

rule GRCh38_bwa_index:
    input:
        "../../reference/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    output:
        touch("../../reference/GRCh38_p14_genomic.ref")
    conda: "../envs/bwa_snake.yml"
    threads: 6
    params:
        mem="5g"
    shell:
        """
        bwa index {input}
        """

rule GRCh38_bwa_map:
    input:
        bwa_index_done = "../../reference/GRCh38_p14_genomic.ref",
        ref="../../reference/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
        seq="../lc_rmd_reads/{metag}_lc_rmd.fq"
    output:
        temp("../GRCh38_mapped/{metag}.bam")
    conda:
        "../envs/bwa_snake.yml"
    threads: 12
    params: mem="5g"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.seq} | \
        samtools view -@ {threads} -Sb > {output}
        """

rule GRCh38_samtools_sort:
    input:
        "../GRCh38_mapped/{metag}.bam"
    output:
        temp("../GRCh38_sorted/{metag}.bam")
    conda:
        "../envs/bwa_snake.yml"
    threads: 12
    params: mem="5g"
    shell:
        "samtools sort -@ {threads} {input} -o {output}"

rule GRCh38_samtools_index:
    input:
        "../GRCh38_sorted/{metag}.bam"
    output:
        temp("../GRCh38_sorted/{metag}.bam.bai")
    conda:
        "../envs/bwa_snake.yml"
    threads: 1
    params: mem="5g"
    shell:
        "samtools index {input}"

rule GRCh38_bam_aln_filter:
    input:
        rscript="../../../scripts/pcoa_sag_abund/frag_recruit_bam_aln_filter.r",
        bam="../GRCh38_sorted/{metag}.bam",
        bai="../GRCh38_sorted/{metag}.bam.bai"
    output:
        mapped=temp("../GRCh38_filtered_mapped/{metag}.bam"),
        unmapped=temp("../GRCh38_filtered_unmapped/{metag}.bam")
    threads: 12
    params:
        mem="5g",
        aln_ani="0.99",
        aln_cov="0.9",
        aln_len="100"
    shell:
        """
        {input.rscript} {input.bam} \
        {params.aln_ani} {params.aln_cov} {params.aln_len} \
        {output.mapped} {output.unmapped}
        """

rule GRCh38_unmapped_bam_sort_name:
    input:
        "../GRCh38_filtered_unmapped/{metag}.bam"
    output:
        temp("../GRCh38_unmapped_sort/{metag}.bam")
    threads: 12
    conda: "../envs/bwa_snake.yml"
    params: mem="5g"
    shell:
        """
        samtools sort -n -@ {threads} {input} -o {output}
        """

rule GRCh38_mapped_bam_sort_name:
    input:
        "../GRCh38_filtered_mapped/{metag}.bam"
    output:
        temp("../GRCh38_mapped_sort/{metag}.bam")
    threads: 12
    conda: "../envs/bwa_snake.yml"
    params: mem="5g"
    shell:
        """
        samtools sort -n -@ {threads} {input} -o {output}
        """

rule GRCh38_unmapped_bam2fq:
    input:
        "../GRCh38_unmapped_sort/{metag}.bam"
    output:
        temp("../GRCh38_unmapped_seq/{metag}.fq")
    threads: 2
    conda: "../envs/bwa_snake.yml"
    params: mem="5g"
    shell:
        """
        samtools fastq -@ {threads} -f 4 {input} > {output}
        """

rule GRCh38_mapped_bam2fq:
    input:
        "../GRCh38_mapped_sort/{metag}.bam"
    output:
        temp("../GRCh38_mapped_seq/{metag}.fq")
    threads: 2
    conda: "../envs/bwa_snake.yml"
    params: mem="5g"
    shell:
        """
        samtools fastq -@ {threads} -F 4 {input} > {output}
        """

rule GRCh38_mapped_fq2fa:
    input:
        "../GRCh38_mapped_seq/{metag}.fq"
    output:
        "../GRCh38_mapped_seq/{metag}.fa"
    threads: 1
    params: mem="1g"
    shell:
        """
        sed -n '1~4s/^@/>/p;2~4p' {input} > {output}
        """

rule bwa_index_ref_genomes:
    input:
        "../../reference/{ref_genomes}.fa.gz"
    output:
        done=touch("../../reference/{ref_genomes}.ref")
    conda: "../envs/bwa_snake.yml"
    threads: 24
    params: mem="20g"
    shell:
        """
        bwa index {input}
        """

rule metag_sag_bwa_map:
    input:
        bwa_index_done = "../../reference/{ref_genomes}.ref",
        ref="../../reference/{ref_genomes}.fa.gz",
        seq="../GRCh38_unmapped_seq/{metag}.fq",
    output:
        temp("../{ref_genomes}/mapped/{metag}.bam")
    conda:
        "../envs/bwa_snake.yml"
    threads: 15
    params: mem="200g"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.seq} | \
        samtools view -@ {threads} -Sb > {output}
        """

rule metag_sag_samtools_sort:
    input:
        "../{ref_genomes}/mapped/{metag}.bam"
    output:
        "../{ref_genomes}/sorted/{metag}.bam"
    conda:
        "../envs/bwa_snake.yml"
    threads: 15
    params: mem="5g"
    shell:
        """
        samtools sort -@ {threads} {input} -o {output}
        """

rule dustmask:
    input:
        "../../low_entropy/test/gd_sags/{sag}.fasta"
    output:
        "../../low_entropy/dustmasker/{sag}.bed"
    threads: 1
    conda:
        "../../low_entropy/envs/blast.yml"
    params: mem="5g"
    shell:
        """
        dustmasker -in {input} -out {output} -outfmt acclist
        """

rule get_long_lcr_bed:
    input:
        bed=expand("../../low_entropy/dustmasker/{sag}.bed", sag=sag_id),
        rscript=("../../../scripts/pcoa_sag_abund/frag_recruit_lcr.r")
    output:
        "../../low_entropy/long_lcr.bed"
    threads: 15
    params: mem="5g"
    shell:
        """
        {input.rscript} {threads}
        """

rule barrnap_bac:
    input:
       "../../low_entropy/sags/bac/{bac}.fasta"
    output:
        "../../low_entropy/barrnap/{bac}.gff"
    threads: 1
    conda:
        "../../low_entropy/envs/barrnap.yml"
    params:
        outdir="../../low_entropy/barrnap",
        mem="5g"
    shell:
        """
        barrnap --kingdom bac {input} | \
            grep -v "##" > {output} ||
        touch ../barrnap/{wildcards.bac}_no_rrna
        """

rule barrnap_arc:
    input:
       "../../low_entropy/sags/arc/{arc}.fasta"
    output:
        "../../low_entropy/barrnap/{arc}.gff"
    threads: 1
    conda:
        "../../low_entropy/envs/barrnap.yml"
    params:
        outdir="../../low_entropy/barrnap",
        mem="5g"
    shell:
        """
        barrnap --kingdom arc {input} | \
            grep -v "##" > {output} ||
        touch ../../low_entropy/barrnap/{wildcards.arc}_no_rrna
        """

rule get_rrna_bed:
    input:
        bac=expand("../../low_entropy/barrnap/{bac}.gff", bac=bac_id),
        arc=expand("../../low_entropy/barrnap/{arc}.gff", arc=arc_id),
        rscript="../../../scripts/pcoa_sag_abund/frag_recruit_rrna.r"
    output:
        "../../low_entropy/rrna.bed"
    threads: 15
    params: mem="5g"
    shell:
        """
        {input.rscript} {threads}
        """

rule samtools_wo_lcr:
    input:
       "../{ref_genomes}/sorted/{metag}.bam",
       "../../low_entropy/long_lcr.bed"
    output:
        "../{ref_genomes}/low_entropy/wo_lcr/{metag}_w_Regions.bam",
        "../{ref_genomes}/low_entropy/wo_lcr/{metag}_wo_Regions.bam"
    threads: 1
    conda:
        "../../low_entropy/envs/samtools.yml"
    params: mem="5g"
    shell:
        """
        samtools view {input[0]} -b -h -o {output[0]} \
            -U {output[1]} -L {input[1]}
        """

rule samtools_wo_lcr_rrna:
    input:
        "../{ref_genomes}/low_entropy/wo_lcr/{metag}_wo_Regions.bam",
        "../../low_entropy/rrna.bed"
    output:
        "../{ref_genomes}/low_entropy/wo_lcr_rrna/{metag}_w_Regions.bam",
        "../{ref_genomes}/low_entropy/wo_lcr_rrna/{metag}_wo_Regions.bam"
    threads: 1
    conda:
        "../../low_entropy/envs/samtools.yml"
    params: mem="5g"
    shell:
        """
        samtools view {input[0]} -b -h -o {output[0]} \
            -U {output[1]} -L {input[1]}
        """

rule samtools_wo_lcr_rrna_index:
    input:
        "../{ref_genomes}/low_entropy/wo_lcr_rrna/{metag}_wo_Regions.bam"
    output:
        "../{ref_genomes}/low_entropy/wo_lcr_rrna/{metag}_wo_Regions.bam.bai"
    conda:
        "../envs/bwa_snake.yml"
    threads: 1
    params: mem="5g"
    shell:
        "samtools index {input}"

rule bam_aln_filter_count:
    input:
        rscript="../../../scripts/pcoa_sag_abund/frag_recruit_bam_aln_filter_count.r",
        bam="../{ref_genomes}/low_entropy/wo_lcr_rrna/{metag}_wo_Regions.bam",
        bai="../{ref_genomes}/low_entropy/wo_lcr_rrna/{metag}_wo_Regions.bam.bai"
    output:
        stat=temp("../{ref_genomes}/stats/{metag}_mapped_filtered.tmp"),
        bam="../{ref_genomes}/aln_filter_bam/{metag}.bam"
    threads: 12
    params:
        mem="5g",
        aln_ani="0.95",
        aln_cov="0",
        aln_len="100"
    shell:
        """
        {input.rscript} {input.bam} \
        {params.aln_ani} {params.aln_cov} {params.aln_len} \
        {output.stat} {output.bam}
        """

# exclude unmapped and secondary alignments
rule count_total:
    input:
        "../{ref_genomes}/low_entropy/wo_lcr_rrna/{metag}_wo_Regions.bam",
        "../{ref_genomes}/low_entropy/wo_lcr_rrna/{metag}_wo_Regions.bam.bai"
    output:
       temp("../{ref_genomes}/stats/{metag}_total_read.tmp")
    conda: "../envs/bwa_snake.yml"
    threads: 1
    params: mem="5g"
    shell:
        """
        samtools view -c {input[0]} > {output}
        """

rule format:
    input:
        total="../{ref_genomes}/stats/{metag}_total_read.tmp",
        mapped="../{ref_genomes}/stats/{metag}_mapped_filtered.tmp"
    output:
        total=temp("../{ref_genomes}/stats/{metag}_total_read.tmp2"),
        mapped=temp("../{ref_genomes}/stats/{metag}_mapped_filtered.tmp2")
    threads: 1
    params: mem="1g"
    shell:
        """
        awk '{{ print FILENAME","$0 }}' {input.total} | \
            sed 's/.*\///g' | sed 's/_.*,/,/' > {output.total}
        sed 's/.*\///g' {input.mapped} > {output.mapped}
        """

# use double bracket for wildcards not subjected to expand function
rule combine:
    input:
        total=expand("../{{ref_genomes}}/stats/{metag}_total_read.tmp2",
            metag=metag_ids),
        mapped=expand("../{{ref_genomes}}/stats/{metag}_mapped_filtered.tmp2",
            metag=metag_ids)
    output:
        total_final="../{ref_genomes}/stats/n_total_read.csv",
        mapped_final="../{ref_genomes}/stats/n_mapped_read.csv"
    threads: 1
    params: mem="1g"
    shell:
        """
        echo "run_accessions,n_total_reads" > {output.total_final}; \
            cat {input.total} >> {output.total_final}
        echo "run_accessions,n_mapped_reads" > {output.mapped_final}; \
            cat {input.mapped} >> {output.mapped_final}
        """

#! below: summarizing the number of post-qc read
rule count_prior_rarefy:
    input:
        "../merged/{metag}.extendedFrags.fastq.gzip"
    output:
        temp("../{ref_genomes}/post_qc_summary/{metag}_read_count.csv")
    params:
        tmpdir="../{ref_genomes}/post_qc_tmp",
        mem="5g"
    threads: 1
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input} | wc -l > {params.tmpdir}/{wildcards.metag}_read_count.txt && \
        awk '{{ print FILENAME","$0/4 }}' {params.tmpdir}/{wildcards.metag}_read_count.txt | \
        sed 's/.*\///g' | sed 's/_.*,/,/' > {output}  && \
        rm {params.tmpdir}/{wildcards.metag}_read_count.txt
        """

rule combine_prior_rarefy:
    input:
        expand("../{{ref_genomes}}/post_qc_summary/{metag}_read_count.csv",
            metag=metag_ids)
    output:
        "../{ref_genomes}/post_qc_summary/metag_post_qc_read.csv"
    threads: 1
    params: mem="1g"
    shell:
        """
        echo "run_accessions,survived_reads" > {output}; \
        cat {input} >> {output}
        """
