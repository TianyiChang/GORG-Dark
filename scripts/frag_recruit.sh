#!/bin/bash

maindir=/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit
codedir=/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/scripts/pcoa_sag_abund
old_codedir=/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/scripts/.vscode

mkdir -p $maindir/result_4_sra_metag/snakefiles
mkdir -p $maindir/result_4_sra_metag/envs
mkdir -p $maindir/result_4_local_metag/snakefiles
mkdir -p $maindir/result_4_local_metag/envs
mkdir -p $maindir/result_4_local_maria_cariaco_metag/snakefiles
mkdir -p $maindir/result_4_local_maria_cariaco_metag/envs
mkdir -p $maindir/result_4_sra_dark_sunlit/snakefiles
mkdir -p $maindir/result_4_sra_dark_sunlit/envs
mkdir -p $maindir/result_4_local_dark_sunlit/snakefiles
mkdir -p $maindir/result_4_local_dark_sunlit/envs
mkdir -p $maindir/result_4_sra_dark/snakefiles
mkdir -p $maindir/result_4_sra_dark/envs
mkdir -p $maindir/result_4_local_dark/snakefiles
mkdir -p $maindir/result_4_local_dark/envs
mkdir -p $maindir/result_4_sra_dark_gdv4omd/snakefiles
mkdir -p $maindir/result_4_sra_dark_gdv4omd/envs

mkdir -p $maindir/metag/from_collab
mkdir -p $maindir/reference
mkdir -p $maindir/metag_wildcards/from_sra
mkdir -p $maindir/metag_wildcards/from_collab

################################
## Update metag or reference ##
################################

# # first manually update "metag_full_list" if adding more metagenomes
# $maindir/metadata/metag_full_list.txt
# #20230903: include additional seas metag (e.g., Black Sea) from "all_dark_sra_run_v2.txt" produced by "metag_metat_metadata.r"

# # then use the code to update wildcards
# for i in $(cat "${maindir}/metadata/metag_full_list.txt"); do
#     filepath="${maindir}/metag_wildcards/${i}"
    
#     if [ ! -e "$filepath" ]; then
#         touch "$filepath"
#     else
#         echo "File '$filepath' already exists and will not be replaced."
#     fi
# done

#! 1. update the below tables without changing file names
head $maindir/metadata/local_metag_list.txt
head $maindir/metadata/sra_run_list.txt

# notes: 20230904: remove SRR4028169 from "sra_run_list.txt", failed to split into PE files
# notes: 20230908: include 12 additional Black Sea metag after adjusting aphotic zone depth, see "metag_metat_metadata.r"
# notes: 20231223: include 66 representative epi metag, see "frag_recruit_select_epi_metag.r".
#todo: After smk run from the above update [20231223]: re-run "frag_recruit_stats_depth_layer_indiv_sags.r",
#todo: then re-run "frag_recruit_pcoa.r", then "generate_pcoa_plots_sag_abund.r"

#! There are only 2 dark ocean metags for Mediterranean Sea
# notes: 20240401: [before] sra_run_list.txt contains 320 dark ocean metag. [new] add 263 dark ocean metag from 2024 malaspina metag, add 278 latitude longitude derep surface metag 
# notes: 20240401: [before] sra_run_list.txt contains 796 nr dark and surface sra run
# notes: 20240626: [latest] sra_run_list_dark_only.txt contains 516 nr dark sra run. [new] add 9 (5 from the hadal zone) dark ocean mariana trench metag, see 'liu_jiwen_2023_microbiome_tables1.xlsx'
# notes: 20240627: [latest] sra_run_list.txt contains 796 nr dark and surface sra run. [new] add 9 (5 from the hadal zone) dark ocean mariana trench metag, see 'liu_jiwen_2023_microbiome_tables1.xlsx'

#! 2. update "ref_genomes" in "frag_recruit.smk" if adding more references
#! add the assemblies into reference folder
# contig names have been unified between Black Sea and others
cat $maindir/../dark_v3/*.fasta | \
    sed 's/^>SCGC_/>/' | \
    sed 's/_contig/_NODE_/' > $maindir/reference/gorg_v3_concat.fa

gzip $maindir/reference/gorg_v3_concat.fa

# 20240401: combine gorg-tropics and gorg-dark
cat $maindir/reference/gorg_v3_concat.fa.gz \
    $maindir/reference/outside_gorg_tropics.fa.gz > $maindir/reference/gorg_v3_tropics.fa.gz

# 20240625: dark_v4 against dark metag only with 263 addit 024 malaspina metag
cat $maindir/../dark_v4/*.fasta | \
    sed 's/^>SCGC_/>/' | \
    sed 's/_contig/_NODE_/' > $maindir/reference/gorg_v4_concat.fa

conda activate pigz
pigz $maindir/reference/gorg_v4_concat.fa -p 12

# 20240627: combine dark_v4 and omd-m
cat $maindir/reference/gorg_v4_concat.fa.gz \
    $maindir/reference/outside_omd_mags.fa.gz > $maindir/reference/gorg_v4_omd.fa.gz

# 20240627: combine gorg-tropics and gorg-dark_v4
cat $maindir/reference/gorg_v4_concat.fa.gz \
    $maindir/reference/outside_gorg_tropics.fa.gz > $maindir/reference/gorg_v4_tropics.fa.gz

# 20240715: combine dark_v4 and omd-m-dark
cat $maindir/reference/gorg_v4_concat.fa.gz \
    $maindir/reference/omd_m_dark.fa.gz > $maindir/reference/gorg_v4_omd_dark.fa.gz


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
conda env export --name sra-tools --file $maindir/result_4_sra_dark_gdv4omd/envs/sra_tools.yml
conda env export --name pigz --file $maindir/result_4_sra_dark_gdv4omd/envs/pigz.yml
#conda env export --name trimmomatic --file $maindir/result_4_sra_dark_gdv4omd/envs/trimmomatic.yml
conda env export --name flash --file $maindir/result_4_sra_dark_gdv4omd/envs/flash.yml
conda env export --name seqtk --file $maindir/result_4_sra_dark_gdv4omd/envs/seqtk.yml
conda env export --name bwa_snake --file $maindir/result_4_sra_dark_gdv4omd/envs/bwa_snake.yml
conda env export --name blast_2.13 --file $maindir/result_4_sra_dark_gdv4omd/envs/blast.yml
conda env export --name barrnap --file $maindir/result_4_sra_dark_gdv4omd/envs/barrnap.yml

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
#! noted 20240401: next time run for the below codes, the metag analyzed will be changed (now it includes malaspina and sunlit ocean metag)
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

#! 3. RUN 4 local maria's Cariaco metagenomes
conda activate snakemake

cp $codedir/frag_recruit_local_maria_cariaco.smk $maindir/result_4_local_maria_cariaco_metag/snakefiles
cd $maindir/result_4_local_maria_cariaco_metag/snakefiles

snakemake -s frag_recruit_local_maria_cariaco.smk --rulegraph | dot -Tpdf > frag_recruit_local_maria_cariaco_dag.pdf

# dry run
snakemake -n -s frag_recruit_local_maria_cariaco.smk --rerun-triggers mtime

# using cluster
# option: --keep-going: do not stop if some of the jobs failed
# option: --rerun-triggers mtime: do not want to re-run jobs after changing contents in some rules)
# option:  --rerun-incomplete

snakemake -s frag_recruit_local_maria_cariaco.smk --use-conda \
--cluster 'qsub -q low -l ncpus={threads},mem={resources.mem_mb}mb,walltime=48:10:00' \
-j 200 --latency-wait 120 --keep-going --rerun-triggers mtime


#! 4. RUN 4 dark and sunlit sra metagenomes with gorgd and gorgt
conda activate snakemake

cp $codedir/frag_recruit_sra_dark_sunlit.smk $maindir/result_4_sra_dark_sunlit/snakefiles
cd $maindir/result_4_sra_dark_sunlit/snakefiles

snakemake -s frag_recruit_sra_dark_sunlit.smk --rulegraph | dot -Tpdf > frag_recruit_sra_dark_sunlit_dag.pdf

# dry run
snakemake -n -s frag_recruit_sra_dark_sunlit.smk --rerun-triggers mtime

# using cluster
# option: --keep-going: do not stop if some of the jobs failed
# option: --rerun-triggers mtime: do not want to re-run jobs after changing contents in some rules)
# option:  --rerun-incomplete

snakemake -s frag_recruit_sra_dark_sunlit.smk --use-conda \
    --cluster 'qsub -q low -l ncpus={threads},mem={resources.mem_mb}mb,walltime=48:10:00' \
    -j 200 --latency-wait 120 --keep-going --rerun-triggers mtime --scheduler greedy


#! 5. RUN 4 dark and sunlit local metagenomes with gorgd and gorgt
conda activate snakemake

cp $codedir/frag_recruit_local_dark_sunlit.smk $maindir/result_4_local_dark_sunlit/snakefiles
cd $maindir/result_4_local_dark_sunlit/snakefiles

snakemake -s frag_recruit_local_dark_sunlit.smk --rulegraph | dot -Tpdf > frag_recruit_sra_localk_sunlit_dag.pdf

# dry run
snakemake -n -s frag_recruit_local_dark_sunlit.smk --rerun-triggers mtime

# using cluster
# option: --keep-going: do not stop if some of the jobs failed
# option: --rerun-triggers mtime: do not want to re-run jobs after changing contents in some rules)
# option:  --rerun-incomplete

snakemake -s frag_recruit_local_dark_sunlit.smk --use-conda \
    --cluster 'qsub -q low -l ncpus={threads},mem={resources.mem_mb}mb,walltime=48:10:00' \
    -j 200 --latency-wait 120 --keep-going --rerun-triggers mtime

#!/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/reference/gorg_v4_omd_dark.ref
#! 6. RUN 4 dark sra metagenomes with dark_v4
# exclude epi metag from sra_run_list.txt
cd $maindir/metadata
grep -v -f selected_epi_metag_runs.txt sra_run_list.txt > sra_run_list_dark_only.txt

conda activate snakemake

cp $codedir/frag_recruit_sra_dark.smk $maindir/result_4_sra_dark/snakefiles
cd $maindir/result_4_sra_dark/snakefiles

# dry run
snakemake -n -s frag_recruit_sra_dark.smk --rerun-triggers mtime

# using cluster
# option: --keep-going: do not stop if some of the jobs failed
# option: --rerun-triggers mtime: do not want to re-run jobs after changing contents in some rules)
# option:  --rerun-incomplete
# option:  --scheduler greedy if there are too many jobs

snakemake -s frag_recruit_sra_dark.smk --use-conda \
    --cluster 'qsub -q low -l ncpus={threads},mem={resources.mem_mb}mb,walltime=48:10:00' \
    -j 200 --latency-wait 120 --keep-going --rerun-triggers mtime --scheduler greedy

#! 7. RUN 4 dark local metagenomes
conda activate snakemake

cp $codedir/frag_recruit_local_dark.smk $maindir/result_4_local_dark/snakefiles
cd $maindir/result_4_local_dark/snakefiles

# dry run
snakemake -n -s frag_recruit_local_dark.smk --rerun-triggers mtime

# using cluster
# option: --keep-going: do not stop if some of the jobs failed
# option: --rerun-triggers mtime: do not want to re-run jobs after changing contents in some rules)
# option:  --rerun-incomplete

snakemake -s frag_recruit_local_dark.smk --use-conda \
    --cluster 'qsub -q low -l ncpus={threads},mem={resources.mem_mb}mb,walltime=48:10:00' \
    -j 200 --latency-wait 120 --keep-going --rerun-triggers mtime

#todo: rerun 'frag_recruit_mapping_rate_across_studies.r'