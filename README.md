# Usage
The main Python script PairedEndVariantConda.py aims at performing variant calling analysis of bacterial genomes from a Python module genomic.py and a Conda environment PairedEndVariantCalling.
- This workflow run reference indexing for BWA/Picard/Samtools (step 1_reference), BBnorn (step 2_normalization), Trimmomatic (step 3_trimming), BWA and Picard (step 4_mapping), GATK4 (step 5_calling) and quality assessment of mapping with Samtools (step 6_quality) based on paired-end reads from a single genomic sample, successively.
- The main script PairedEndVariantConda.py and module genomic.py (version 20201006, October 2020) were prepared and tested with Python and dependencies below.
- The module genomic.py has to be with the present main script PairedEndVariantConda.py to launch it properly.
- The Conda environment PairedEndVariantCalling has to be prepared as presented below.
- The user can setup his own dependencies in his own bin.
- The paired-end reads must be named ID_R1.fastq.gz and ID_R2.fastq.gz for forward and reverse reads, respectively (ID means sample identifier).
- The IDs have to include a maximum of 16 alphanumeric characters (AZ, az , 09) including potential underscores (_).
- The accents (‘, ¨, ^), space ( ), hyphen (-), and special characters (/, '\', », (, }, =, +, @) are not accepted in the IDs.
- The quality scores paired-end reads must be encoded with Phred33.
# Dependencies
The main script PairedEndVariantConda.py and module genomic.py (version 20201006) were prepared and tested with Conda packages below (Name/Version/Build/Channel).
- bwa/0.7.17/hed695b0_7/bioconda
- samtools/1.10/h2e538c0_3/bioconda
- picard/2.23.4/0/bioconda
- bbmap/38.86/h1296035_0/bioconda
- trimmomatic/0.39/1/bioconda
- python/3.8.5/h1103e12_9_cpython/conda-forge
- biopython/1.78/py38h1e0a361_0/conda-forge
- gatk4/4.1.8.1/py38_0/bioconda
# Building of the Conda Environment PairedEndVariantCalling
## 1/ From available targeted Conda packages
```
conda activate
conda create -n PairedEndVariantCalling
conda activate PairedEndVariantCalling
conda search bwa
conda install -c bioconda bwa=0.7.17=hed695b0_7
conda search samtools
conda install -c bioconda samtools=1.10=h2e538c0_3
conda search picard
conda install -c bioconda picard=2.23.4=0
conda search bbmap
conda install -c bioconda bbmap=38.86=h1296035_0
conda search trimmomatic
conda install -c bioconda trimmomatic=0.39=1
conda search python
conda install -c conda-forge python=3.8.5=h1103e12_9_cpython
conda search biopython
conda install -c conda-forge biopython=1.78=py38h1e0a361_0
conda search gatk4
conda install -c bioconda gatk4=4.1.8.1=py38_0
```
## 2/ From available updated Conda packages
```
conda activate
conda create -n PairedEndVariantCalling
conda activate PairedEndVariantCalling
conda install -c bioconda bwa
conda update -c bioconda bwa
conda install -c bioconda samtools
conda update -c bioconda samtools
conda install -c bioconda picard
conda update -c bioconda picard
conda install -c bioconda bbmap
conda update -c bioconda bbmap
conda install -c bioconda trimmomatic
conda update -c bioconda trimmomatic
conda install -c conda-forge python
conda update -c conda-forge python
conda install -c conda-forge biopython
conda update -c conda-forge biopython
conda install -c bioconda gatk4
conda update -c bioconda gatk4
```
# Launching of the script PairedEndVariantConda.py
## 1/ With a single set of paired-end reads
### 1.1/ prepare a single command in a Bash script (bash_PairedEndVariantConda.sh)
```
#!/bin/bash
#SBATCH -p Research
#SBATCH -o %x.%N.%j.out
#SBATCH -e %x.%N.%j.err
#SBATCH --cpus-per-task=48
#SBATCH --job-name=test-20201005
source /global/conda/bin/activate;conda activate PairedEndVariantCalling; \
python /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/PairedEndVariantConda.py \
 -r VariantCalling \
 -t 40 \
 -n 100 \
 -R1 /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/ERR3997409_R1.fastq.gz \
 -R2 /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/ERR3997409_R2.fastq.gz \
 -adap /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/NexteraPE-PE.fa \
 -ref /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/Enteritidis_P125109/Enteritidis_P125109.fasta \
 -norm /global/conda/envs/PairedEndVariantCalling/bin/bbnorm.sh \
 -trim /global/conda/envs/PairedEndVariantCalling/bin/trimmomatic \
 -align /global/conda/envs/PairedEndVariantCalling/bin/bwa \
 -process /global/conda/envs/PairedEndVariantCalling/bin/picard \
 -call /global/conda/envs/PairedEndVariantCalling/bin/gatk \
 -cov /global/conda/envs/PairedEndVariantCalling/bin/samtools
```
### 1.2/ run the Bash script bash_PairedEndVariantConda.sh with sbatch
```
sbatch bash_PairedEndVariantConda.sh
```
### 1.3/ compile the depth and breadth coverages from mapping
```
grep . /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/6_quality/*.bam.* > /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/6_quality/coverage.metrics
```
## 2/ With multiple sets of paired-end reads
### 2.1/ creat a file list_of_IDs.lst including a list of ID samples to process (one ID per line with \n)
```
gedit list_of_IDs.lst
```
### 2.2/ creat a file commands.lst including a list of Bash commands
```
rm commands.lst
for l in `cat list_of_IDs.lst`; do
	echo "source /global/conda/bin/activate;conda activate PairedEndVariantCalling;python /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/PairedEndVariantConda.py -r VariantCalling -t 40 -n 100 -R1 /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/$l _R1.fastq.gz -R2 /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/$l _R2.fastq.gz -adap /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/NexteraPE-PE.fa -ref /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/Enteritidis_P125109/Enteritidis_P125109.fasta -norm /global/conda/envs/PairedEndVariantCalling/bin/bbnorm.sh -trim /global/conda/envs/PairedEndVariantCalling/bin/trimmomatic -align /global/conda/envs/PairedEndVariantCalling/bin/bwa -process /global/conda/envs/PairedEndVariantCalling/bin/picard -call /global/conda/envs/PairedEndVariantCalling/bin/gatk -cov /global/conda/envs/PairedEndVariantCalling/bin/samtools";
done >> commands.lst
sed -i "s@ _R1.fastq.gz@_R1.fastq.gz@" commands.lst
sed -i "s@ _R2.fastq.gz@_R2.fastq.gz@" commands.lst
```
### 2.3/ run the Bash commands of the file commands.lst with sarray
```
sarray -p Research --cpus-per-task=48 -e %x.%N.%j.err -o %x.%N.%j.out --job-name=test-20201006 commands.lst
```
### 2.4/ compile the depth and breadth coverages from mapping
```
grep . /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/6_quality/*.bam.* > /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/6_quality/coverage.metrics
```
# Acknowledgment
My old colleagues Arnaud Felten and Ludovic Mallet with whom I learned a lot
# Author
Nicolas Radomski
