#!/bin/bash

maindir=/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit
codedir=/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/scripts/pcoa_sag_abund
old_codedir=/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/scripts/.vscode

mkdir -p $maindir/result_4_sra_metag/snakefiles
mkdir -p $maindir/result_4_sra_metag/envs
mkdir -p $maindir/result_4_local_metag/snakefiles
mkdir -p $maindir/result_4_local_metag/envs

mkdir -p $maindir/metag/from_collab
mkdir -p $maindir/reference
mkdir -p $maindir/metag_wildcards/from_sra
mkdir -p $maindir/metag_wildcards/from_collab

################################
## Update metag or reference ##
################################

# # first manually update "metag_full_list" if adding more metagenomes
# $maindir/metadata/metag_full_list.txt
# #20230903: include additional seas metag (e.g., Black Sea) from "all_dark_sra_run_v2.txt" produced by "metadata.r"

# # then use the code to update wildcards
# for i in $(cat "${maindir}/metadata/metag_full_list.txt"); do
#     filepath="${maindir}/metag_wildcards/${i}"
    
#     if [ ! -e "$filepath" ]; then
#         touch "$filepath"
#     else
#         echo "File '$filepath' already exists and will not be replaced."
#     fi
# done

#! 1. manually modify the below tables without changing file names
head $maindir/metadata/local_metag_list.txt
head $maindir/metadata/sra_run_list.txt

notes: 20230904: remove SRR4028169 from "sra_run_list.txt", failed to split into PE files
notes: 20230908: include 12 additional Black Sea metag after adjusting aphotic zone depth, see "metadata.r"

#! 2. update "ref_genomes" in "frag_recruit.smk" if adding more references
#! add the assemblies into reference folder
# contig names have been unified between Black Sea and others
cat $maindir/../dark_v3/*.fasta | \
    sed 's/^>SCGC_/>/' | \
    sed 's/_contig/_NODE_/' > $maindir/reference/gorg_v3_concat.fa

gzip $maindir/reference/gorg_v3_concat.fa

########################################
# re-organizing file and code systems ##
########################################

cp $maindir/mapping/addit_collaboraters/fastq/*.fq.gz \
    $maindir/metag/from_collab

(cp $maindir/mapping/reference/* $maindir/reference) && rm $maindir/reference/gorg_v2_*

cp $old_codedir/frag_recruit_bam_aln_filter_count.r $codedir
cp $old_codedir/frag_recruit_bam_aln_filter.r $codedir
cp $old_codedir/frag_recruit_rrna.r $codedir
cp $old_codedir/frag_recruit_bam_aln_filter_count.r $codedir
cp $old_codedir/frag_recruit_lcr_reads.r $codedir

for i in $(ls $maindir/metag/from_collab/*_R1.fq.gz | sed 's/_R1.fq.gz//'); do
    mv ${i}_R1.fq.gz ${i}_1.fastq.gz;
done

for i in $(ls $maindir/metag/from_collab/*_R2.fq.gz | sed 's/_R2.fq.gz//'); do
    mv ${i}_R2.fq.gz ${i}_2.fastq.gz;
done

# create conda envs
conda env export --name sra-tools --file $maindir/result_4_metag/envs/sra_tools.yml
conda env export --name pigz --file $maindir/result_4_metag/envs/pigz.yml
#conda env export --name trimmomatic --file $maindir/result_4_metag/envs/trimmomatic.yml
conda env export --name flash --file $maindir/result_4_metag/envs/flash.yml
conda env export --name seqtk --file $maindir/result_4_metag/envs/seqtk.yml
conda env export --name bwa_snake --file $maindir/result_4_metag/envs/bwa_snake.yml
conda env export --name blast_2.13 --file $maindir/result_4_metag/envs/blast.yml
conda env export --name barrnap --file $maindir/result_4_metag/envs/barrnap.yml

###############
## SNAKEMAKE ##
###############

#! 1. RUN 4 local metagenomes
conda activate snakemake

cp $codedir/frag_recruit_local.smk $maindir/result_4_local_metag/snakefiles
cd $maindir/result_4_local_metag/snakefiles

snakemake -s frag_recruit_local.smk --rulegraph | dot -Tpdf > frag_recruit_local_dag.pdf

# dry run
snakemake -n -s frag_recruit_local.smk

# using cluster
# option: --keep-going: do not stop if some of the jobs failed
# option: --rerun-triggers mtime: do not want to re-run jobs after changing contents in some rules)
# option:  --rerun-incomplete

snakemake -s frag_recruit_local.smk --use-conda \
--cluster 'qsub -q low -l ncpus={threads},mem={params.mem},walltime=48:10:00' \
-j 200 --latency-wait 120 --keep-going

#! 2. RUN 4 sra metagenomes
conda activate snakemake

cp $codedir/frag_recruit_sra.smk $maindir/result_4_sra_metag/snakefiles
cd $maindir/result_4_sra_metag/snakefiles

snakemake -s frag_recruit_sra.smk --rulegraph | dot -Tpdf > frag_recruit_sra_dag.pdf

# dry run
snakemake -n -s frag_recruit_sra.smk --rerun-triggers mtime

# using cluster
# option: --keep-going: do not stop if some of the jobs failed
# option: --rerun-triggers mtime: do not want to re-run jobs after changing contents in some rules)
# option:  --rerun-incomplete

snakemake -s frag_recruit_sra.smk --use-conda \
--cluster 'qsub -q low -l ncpus={threads},mem={resources.mem_mb},walltime=48:10:00' \
-j 200 --latency-wait 120 --keep-going --rerun-triggers mtime
