Build 1000G reference panel
===========================

Snakemake pipeline for building 1000G reference panel. No filtering is done apart from removing one pair of related participants.

Possible filters to consider:
* Remove low frequency markers (e.g., <5 counts)
* Remove structural variants
* Remove non-unique identifiers


Overview
--------

1. Download the 1000G phase 3 data.

2. Update rsids with latest dbSNP (rsids missing from sex and MT chromosomes).

For each super population of the 1000G phase 3 data on both hg19 and hg38 coordinate systems:

3. Create binary plink files per chromosome.

4. [Optional] Create concatenated autosome vcf and plink file.


Quickstart
-----------

```bash
# clone the repo
git clone https://github.com/letaylor/reference_1000G

# set code base path
SNK_REPO="`pwd`/reference_1000G"
```

If you are running jobs on the cluster, it is best to first start a [tmux](https://github.com/tmux/tmux) session so that the session may be re-attached at a later time point.

```bash
# start session
tmux new -s mk_1000G

# log back into a session
tmux -a mk_1000G
```

Download the required data.

```bash
# make a directory for a demo
mkdir 1000G_phase3
cd 1000G_phase3

# set up the directory with the required reference_data dir
mkdir reference_data

# Download the latest dbSNP information
# see here https://www.ncbi.nlm.nih.gov/dbvar/content/org_summary/
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi
mv All_20180423.vcf.gz reference_data/hg19-dbsnp.vcf.gz
mv All_20180423.vcf.gz.tbi reference_data/hg19-dbsnp.vcf.gz.tbi

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi
mv All_20180418.vcf.gz reference_data/hg38-dbsnp.vcf.gz
mv All_20180418.vcf.gz.tbi reference_data/hg38-dbsnp.vcf.gz.tbi


# Download the 1000 G data
wget --no-parent -r ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
DAT_1000G="$(pwd)/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"

# get sample information (sex and population)
cp ${DAT_1000G}/integrated_call_samples_v3.20130502.ALL.panel reference_data/sample_info.tsv
cp ${DAT_1000G}/20140625_related_individuals.txt reference_data/sample_info-related.tsv

# set up unprocessed reference_data dir
# ... for hg19
mkdir reference_data/hg19-all_samples
for i in $(ls ${DAT_1000G}/*vcf.gz); do
    new_file_name=$(basename $i | sed s/".phase3"/"\t"/g | cut -f 1 | sed s/"ALL."//)
    ln -s $i data/hg19-all_samples/unprocessed-$new_file_name.vcf.gz
done

for i in $(ls ${DAT_1000G}/*vcf.gz.tbi); do
    new_file_name=$(basename $i | sed s/".phase3"/"\t"/g | cut -f 1 | sed s/"ALL."//)
    ln -s $i reference_data/hg19-all_samples/unprocessed-$new_file_name.vcf.gz.tbi
done

# ... for hg38
mkdir reference_data/hg38-all_samples
for i in $(ls ${DAT_1000G}/supporting/GRCh38_positions/*genotypes*vcf.gz); do
    new_file_name=$(basename $i | sed s/"_GRCh38"/"\t"/g | cut -f 1 | sed s/"ALL."//)
    ln -s $i reference_data/hg38-all_samples/unprocessed-$new_file_name.vcf.gz
done

for i in $(ls ${DAT_1000G}/supporting/GRCh38_positions/*genotypes*vcf.gz.tbi); do
    new_file_name=$(basename $i | sed s/"_GRCh38"/"\t"/g | cut -f 1 | sed s/"ALL."//)
    ln -s $i reference_data/hg38-all_samples/unprocessed-$new_file_name.vcf.gz.tbi
done

# set up the config file inside the working dir
# change json file to fit the analysis
cp ${SNK_REPO}/config_analysis.json .

# now run snakemake in dryrun
snakemake --snakefile ${SNK_REPO}/Snakefile --configfile config_analysis.json --dryrun

# now run snakemake
snakemake --snakefile ${SNK_REPO}/Snakefile --configfile config_analysis.json --printshellcmds
```

The above jobs may take a little while, it would be much faster if we utilized a cluster to run the jobs.

Below is a demo for SLURM, you can easily tweak the config file for other systems.

```bash
# set up the config file inside the working dir
# change json files to fit the system
cp ${SNK_REPO}/config_cluster_slurm.json config_cluster.json

# now run snakemake below is an example for SLURM
snakemake --snakefile ${SNK_REPO}/Snakefile --configfile config_analysis.json --printshellcmds --latency-wait 600 --jobs 999 --cluster-config config_cluster.json --cluster 'sbatch --job-name="{cluster.name}" --export=ALL --ntasks={cluster.tasks} --cpus-per-task={cluster.cpus} --mem={cluster.memory}G --output={cluster.output} --error={cluster.error}'

# alternatively use the cluster script (which will load settings from the config file)
snakemake --snakefile ${SNK_REPO}/Snakefile --configfile config_analysis.json --printshellcmds --latency-wait 600 --jobs 999 --cluster-config config_cluster.json --cluster ${SNK_REPO}/wrappers/cluster/slurm.py
```


Reproducibility: Conda   
----------------------

One option to enhance reproducibility is to install software used via Conda. Note that some software (e.g., LocusZoom) will still be missing since they are not currently in a Conda channel.

```bash
# install environment using Conda
conda env create --name mk_1000G --file ${SNK_REPO}/environment.yml

# activate the new Conda environment
source activate mk_1000G
```
