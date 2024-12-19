maindir=/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/related_studies/omd
codedir=/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/scripts/pcoa_sag_abund
outdir=$maindir/omd_m_dark

# Download and decompress OMD genomes
conda activate aria2

aria2c https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/genomes-fasta.tar.gz

conda activate pigz
pigz -d genomes-fasta.tar.gz | tar xf && rm genomes-fasta.tar.gz

# Remove GORG MARD and non-bacteria non-archaea from the downloaded genomes
chmod +x $codedir/subset_omd_mags.r
$codedir/subset_omd_mags.r

mkdir -p $outdir/omd_mags
mkdir -p $outdir/metad/snakefiles

for mag in $(cat mag_list.txt); do
    cp genomes-fasta/$mag omd_m_dark/omd_mags
done

# Retrieve metadata of BioSample where the MAGs recovered
cd $outdir/metad/snakefiles

conda activate snakemake

# dry run
snakemake -n -s $codedir/get_dark_ocean_omd_mags.smk

snakemake -s $codedir/get_dark_ocean_omd_mags.smk --cores 30 --use-conda --resources load=100 --keep-going

# esearch -db BioSample -query SAMEA2620256[Accession] | efetch -format xml > test.xml

# Get OMD MAGs recovered from dark ocean metag
chmod +x $codedir/get_omd_dark.r

$codedir/get_omd_dark.r

# Concat dark ocean OMD-M
for file in $(cat $maindir/omd_m_dark/omd_dark.txt); do
    cat $maindir/omd_m_dark/omd_mags/$file >> $maindir/../../frag_recruit/reference/omd_m_dark.fa
done

cd $maindir/../../frag_recruit/reference
conda activate pigz
pigz omd_m_dark.fa
