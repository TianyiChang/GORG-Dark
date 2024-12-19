#!/bin/sh

# Retrieve metadata of SRA used in frag recruit analyses

maindir=/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/metadata/metag_metadata_240420
codedir=/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/scripts/pcoa_sag_abund

mkdir -p $maindir/snakefiles


conda activate snakemake

cp $codedir/accession2metadata_ncbi.smk $maindir/snakefiles
cd $maindir/snakefiles

# snakemake -s accession2metadata_ncbi.smk --rulegraph | dot -Tpdf > accession2metadata_ncbi_dag.pdf

# dry run
snakemake -n -s accession2metadata_ncbi.smk --rerun-triggers mtime

# using cluster
# option: --keep-going: do not stop if some of the jobs failed
# option: --rerun-triggers mtime: do not want to re-run jobs after changing contents in some rules)
# option:  --rerun-incomplete

# use load to customize max number of jobs by a rule
snakemake -s accession2metadata_ncbi.smk --cores 30 --use-conda --resources load=100






# Install entrez-direct conda pck
mamba create -n entrez bioconda::entrez-direct

# Retrieve metadata
conda activate entrez

cd $maindir

# for SRA Accession
esearch -db sra -query SRR2103008[Accession] | efetch -format xml | grep 'PRIMARY_ID'

# for BioSample
# esearch -db sra -query 'SAMEA2620824[BioSample]' | efetch -format xml




# # Download the complete NCBI SRA metadata
# wget https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/NCBI_SRA_Metadata_Full_20240321.tar.gz

