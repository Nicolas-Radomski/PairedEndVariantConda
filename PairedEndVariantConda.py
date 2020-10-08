#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

#### author: Nicolas Radomski ####
# version for conda in a cluster: conda enviroment PairedEndVariantCalling
# run reference indexing for BWA/Picard/Samtools (step 1_reference), BBnorn (step 1_normalization), Trimmomatic (step 2_trimming), BWA and Picard (4_mapping), GATK4 (5_calling) and quality assessment of mapping with Samtools (6_quality) based on paired-end reads from a single genomic sample
# the module genomic.py has to be with the present main script PairedEndAssemblyConda.py to lunch it properly
# the present main script PairedEndVariantConda.py and module genomic.py (version 20201006, Octobre 2020) were prepared and tested with Python and Conda packages below (Name/Version/Build/Channel)
#- bwa/0.7.17/hed695b0_7/bioconda
#- samtools/1.10/h2e538c0_3/bioconda
#- picard/2.23.4/0/bioconda
#- bbmap/38.86/h1296035_0/bioconda
#- trimmomatic/0.39/1/bioconda
#- python/3.8.5/h1103e12_9_cpython/conda-forge
#- biopython/1.78/py38h1e0a361_0/conda-forge
#- gatk4/4.1.8.1/py38_0/bioconda
# the present main script PairedEndAssemblyVariant.py executes more precisly the commands below
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/bwa index -a is /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/1_reference/Enteritidis_P125109.fasta
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/samtools faidx /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/1_reference/Enteritidis_P125109.fasta
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/picard CreateSequenceDictionary REFERENCE=/global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/1_reference/Enteritidis_P125109.fasta OUTPUT=/global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/1_reference/Enteritidis_P125109.dict
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/bbnorm.sh in=/global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/ERR3997409_R1.fastq.gz in2=/global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/ERR3997409_R2.fastq.gz out=/global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/2_normalization/ERR3997409_R1_N.fastq.gz out2=/global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/2_normalization/ERR3997409_R2_N.fastq.gz target=100 threads=40
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/trimmomatic PE -threads 40 -phred33 /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/2_normalization/ERR3997409_R1_N.fastq.gz /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/2_normalization/ERR3997409_R2_N.fastq.gz /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/3_trimming/ERR3997409_R1_P.fastq.gz /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/3_trimming/ERR3997409_R1_UP.fastq.gz /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/3_trimming/ERR3997409_R2_P.fastq.gz /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/3_trimming/ERR3997409_R2_UP.fastq.gz ILLUMINACLIP:/global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/bwa mem -t 40 -K 10000000 /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/1_reference/Enteritidis_P125109.fasta /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/3_trimming/ERR3997409_R1_P.fastq.gz /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/3_trimming/ERR3997409_R2_P.fastq.gz > /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/4_mapping/ERR3997409.sam
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/picard SortSam -I /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/4_mapping/ERR3997409.sam -O /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/4_mapping/ERR3997409.sorted.bam --SORT_ORDER coordinate --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/picard AddOrReplaceReadGroups -I /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/4_mapping/ERR3997409.sorted.bam -O /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/4_mapping/ERR3997409.ARRG.bam --SORT_ORDER coordinate --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true --RGID ERR3997409 --RGLB ERR3997409 --RGPU ERR3997409 --RGPL illumina --RGSM ERR3997409
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/picard MarkDuplicates -I /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/4_mapping/ERR3997409.ARRG.bam -O /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/4_mapping/ERR3997409.dedupped.bam -M /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/4_mapping/ERR3997409.dedupped.metrics --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT --REMOVE_SEQUENCING_DUPLICATES true --REMOVE_DUPLICATES true
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/gatk --java-options "-Xms50g -Xmx50g" HaplotypeCaller -R /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/1_reference/Enteritidis_P125109.fasta -I /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/4_mapping/ERR3997409.dedupped.bam -O /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/5_calling/ERR3997409.g.vcf.gz -ERC GVCF -stand-call-conf 30 -ploidy 1 -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/samtools depth /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/4_mapping/ERR3997409.dedupped.bam |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average Depth Coverage (X) =",sum/NR; print "Standard Deviation Depth Coverage (X) =",sqrt(sumsq/NR - (sum/NR)**2)}' > /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/6_quality/ERR3997409.bam.depth
#### Execute /global/conda/envs/PairedEndVariantCalling/bin/samtools coverage /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/4_mapping/ERR3997409.dedupped.bam > /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/6_quality/ERR3997409.bam.metrics
#### Execute cat /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/6_quality/ERR3997409.bam.metrics | awk ' { print $6 } ' | tr -d "\n" | sed 's@coverage@Breadth Coverage (%) = @' > /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/6_quality/ERR3997409.bam.breadth

'''
#### exemple of Bash command (bash_PairedEndVariantConda.sh) ####
#!/bin/bash
#SBATCH -p Research
#SBATCH -o %x.%N.%j.out
#SBATCH -e %x.%N.%j.err
#SBATCH --cpus-per-task=48
#SBATCH --job-name=test-20201006
source /global/conda/bin/activate;conda activate PairedEndVariantCalling; \
python /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/PairedEndVariantConda.py \
 -r VariantCalling \
 -t 48 \
 -n 100 \
 -R1 /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/ERR3997409_R1.fastq.gz \
 -R2 /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/ERR3997409_R2.fastq.gz \
 -adap /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/NexteraPE-PE.fa \
 -ref /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/data/Enteritidis_P125109/Enteritidis_P125109.fasta \
 -norm /global/conda/envs/PairedEndAssembly/bin/bbnorm.sh \
 -trim /global/conda/envs/PairedEndAssembly/bin/trimmomatic \
 -align /global/conda/envs/PairedEndVariantCalling/bin/bwa \
 -process /global/conda/envs/PairedEndVariantCalling/bin/picard \
 -call /global/conda/envs/PairedEndVariantCalling/bin/gatk \
 -cov /global/conda/envs/PairedEndVariantCalling/bin/samtools

#### exemple of Bash command execution ####
sbatch bash_PairedEndVariantConda.sh
'''

import os, sys
import argparse
import genomic

# parse arguments
def get_parser():
	
	# function asking arguments
	parser = argparse.ArgumentParser(description="perform reference indexing (.fasta, .dict, .amb, .ann, .bwt, .fai, .pac, .sa), read normalization (.fastq.gz), read trimming (.fastq.gz), read mapping (.bam), read processing including reads sorting, read group setting and duplication removing (.bam), then variant calling (g.vcf format) and quality assessment of mapping (.bam.breadth and .bam.depth) from forward (ID_R1.fastq.gz) and reverse (ID_R2.fastq.gz) paired-end reads (offset of reads must be Phred 33 and only uppercase or lowercase alphanumeric characters are expected in the ID sample)")

	# setting of arguments

	parser.add_argument('-r', action="store", dest='run',
					type=str, required=True, 
					help='name of the run (REQUIRED)')

	parser.add_argument('-t', action="store", dest='threads',
					type=int, required=True, 
					help='number of threads (REQUIRED)')

	parser.add_argument('-n', action="store", dest='depth',
					type=int, required=True, 
					help='read normalization depth (REQUIRED)')

	parser.add_argument('-R1', action="store", dest='forward',
					type=str, required=True, 
					help='path to the forward read (REQUIRED)')

	parser.add_argument('-R2', action="store", dest='reverse',
					type=str, required=True, 
					help='path to the reverse read (REQUIRED)')

	parser.add_argument('-adap', action="store", dest='adapters',
					type=str, required=True, 
					help='path to adapters for trimming (REQUIRED)')

	parser.add_argument('-ref', action="store", dest='reference',
					type=str, required=True, 
					help='path to reference fasta file (REQUIRED)')

	parser.add_argument('-norm', action="store", dest='bbnorm',
					type=str, required=True, 
					help='path to BBnorm (REQUIRED)')					

	parser.add_argument('-trim', action="store", dest='trimmomatic',
					type=str, required=True, 
					help='path to Trimmomatic (REQUIRED)')

	parser.add_argument('-align', action="store", dest='bwa',
					type=str, required=True, 
					help='path to BWA (REQUIRED)')

	parser.add_argument('-process', action="store", dest='picard',
					type=str, required=True, 
					help='path to Picard (REQUIRED)')

	parser.add_argument('-call', action="store", dest='gatk',
					type=str, required=True, 
					help='path to GATK4 (REQUIRED)')

	parser.add_argument('-cov', action="store", dest='samtools',
					type=str, required=True, 
					help='path to Samtools (REQUIRED)')

	return parser

# ask run, threads, depth of read normalization, R1 path, R2 path, adapter path, 
# BBnorm path, Trimmomatic path, BWA path, Picard path, GATK4 path, Samtools path, 
# then return directories 1_reference, 2_normalization, 3_trimming, 4_mapping, 5_calling and 6_quality

def main():
	
	# get parser object
	parser=get_parser()

	# print parser.help if there are no arguments in the command
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)

	# extract arguments from parser
	Arguments=parser.parse_args()

	# define variables of the arguments	
	r=Arguments.run
	t=Arguments.threads
	n=Arguments.depth
	R1=Arguments.forward
	R2=Arguments.reverse
	AD=Arguments.adapters
	RE=Arguments.reference
	BB=Arguments.bbnorm
	TR=Arguments.trimmomatic
	BW=Arguments.bwa
	PI=Arguments.picard
	GA=Arguments.gatk
	SA=Arguments.samtools

	# get and print working directory (wd)
	wd = os.getcwd()
	print ("#### The current working directory is %s" % wd)

	# extract sample name and ID
	R1name = os.path.basename(R1)
	print("#### The forward read is named %s" %R1name)
	R2name = os.path.basename(R2)
	print("#### The reverse read is named %s" %R2name)
	R1id = R1name.replace('_R1.fastq.gz','')
	print("#### The forward ID is %s" %R1id)
	R2id = R2name.replace('_R2.fastq.gz','')
	print("#### The reverse ID is %s" %R2id)

	# check IDs identity, then prepare output directories and run reference indexing
	if R1id != R2id:
		print("#### The forward (ID_R1.fastq.gz) and reverse (ID_R2.fastq.gz) reads do not match")
		sys.exit("#### Please, check the format of the paired-end reads in {} and R2 in {}".format(R1, R2))
	else:
		print("#### The forward (ID_R1.fastq.gz) and reverse (ID_R2.fastq.gz) reads match")
		# create an output directory called by the run name if it does not exist with the function nodir_makedir_warning of the module genomic.py
		genomic.nodir_makedir_warning(directory = r)
		# create output directories if they do not exist with the function nodir_makedir_warning of the module genomic.py
		genomic.nodir_makedir_warning(directory = r + '/' + '1_reference')
		genomic.nodir_makedir_warning(directory = r + '/' + '2_normalization')
		genomic.nodir_makedir_warning(directory = r + '/' + '3_trimming')
		genomic.nodir_makedir_warning(directory = r + '/' + '4_mapping')
		genomic.nodir_makedir_warning(directory = r + '/' + '5_calling')
		genomic.nodir_makedir_warning(directory = r + '/' + '6_quality')
		# prepare path of output directories
		referenceoutput = wd + '/' + r + '/' + '1_reference' + '/'
		normalizationoutput = wd + '/' + r + '/' + '2_normalization' + '/'
		trimmingoutput = wd + '/' + r + '/' + '3_trimming' + '/'
		mappingoutput = wd + '/' + r + '/' + '4_mapping' + '/'
		callingoutput = wd + '/' + r + '/' + '5_calling' + '/'
		qualityoutput = wd + '/' + r + '/' + '6_quality' + '/'

	# prepare and run reference indexing if the step 1_reference is not done (only of the first processed sample)
	if len(os.listdir(referenceoutput)) > 0:
		print("#### The reference genome is indexed in %s" %referenceoutput)
	else:
		cmdfasta = 'cp' + ' ' + RE + ' ' + referenceoutput
		REname = os.path.basename(RE)
		print("#### The reference is named %s" %REname)
		REid = REname.replace('.fasta','')
		print("#### The reference ID is %s" %REid)
		print("#### Execute %s" %cmdfasta)
		os.system(cmdfasta)
		cmdindexbwa = BW + ' index -a is ' + referenceoutput + REname
		print("#### Execute %s" %cmdindexbwa)
		os.system(cmdindexbwa)
		cmdindexsamtools = SA + ' faidx ' + referenceoutput + REname
		print("#### Execute %s" %cmdindexsamtools)
		os.system(cmdindexsamtools)
		cmdindexpicard = PI + ' CreateSequenceDictionary ' + 'REFERENCE=' + referenceoutput + REname + ' OUTPUT=' + referenceoutput + REid + '.dict'
		print("#### Execute %s" %cmdindexpicard)
		os.system(cmdindexpicard)

	# check the absence of the reference index files of the step 1_reference, then return warning or successfull message with the function absentefile_warning_success of the module genomics.py
	genomic.absentefile_warning_success(expectedfile = r + '/' + '1_reference' + '/' + R1id + '.dict')
	genomic.absentefile_warning_success(expectedfile = r + '/' + '1_reference' + '/' + R1id + '.fasta')
	genomic.absentefile_warning_success(expectedfile = r + '/' + '1_reference' + '/' + R1id + '.fasta.amb')
	genomic.absentefile_warning_success(expectedfile = r + '/' + '1_reference' + '/' + R1id + '.fasta.ann')
	genomic.absentefile_warning_success(expectedfile = r + '/' + '1_reference' + '/' + R1id + '.fasta.bwt')
	genomic.absentefile_warning_success(expectedfile = r + '/' + '1_reference' + '/' + R1id + '.fasta.fai')
	genomic.absentefile_warning_success(expectedfile = r + '/' + '1_reference' + '/' + R1id + '.fasta.pac')
	genomic.absentefile_warning_success(expectedfile = r + '/' + '1_reference' + '/' + R1id + '.fasta.sa')

	# prepare and run BBnorm
	R1N = normalizationoutput + R1id + '_R1_N.fastq.gz'
	R2N = normalizationoutput + R2id + '_R2_N.fastq.gz'
	cmdBBnorm = BB + ' in=' + R1 + ' in2=' + R2 + ' out=' + R1N + ' out2=' + R2N + ' target=' + str(n) + ' threads=' + str(t)
	print("#### Execute %s" %cmdBBnorm)
	os.system(cmdBBnorm)

	# check the absence of the _R1_N.fastq.gz and _R2_N.fastq.gz files of the step 2_normalization, then return warning or successfull message with the function absentefile_warning_success of the module genomics.py
	genomic.absentefile_warning_success(expectedfile = r + '/' + '2_normalization' + '/' + R1id + '_R1_N.fastq.gz')
	genomic.absentefile_warning_success(expectedfile = r + '/' + '2_normalization' + '/' + R1id + '_R2_N.fastq.gz')

	# prepare and run Trimmomatic
	R1P = trimmingoutput + R1id + '_R1_P.fastq.gz'
	R1UP  = trimmingoutput + R1id + '_R1_UP.fastq.gz'
	R2P = trimmingoutput + R2id + '_R2_P.fastq.gz'
	R2UP = trimmingoutput + R2id + '_R2_UP.fastq.gz'
	cmdTrimmomatic = TR + ' PE ' + '-threads ' + str(t) + ' -phred33 ' + R1N + ' ' + R2N + ' ' + R1P + ' ' + R1UP + ' ' + R2P + ' ' + R2UP + ' ILLUMINACLIP:' + AD + ':2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
	print("#### Execute %s" %cmdTrimmomatic)
	os.system(cmdTrimmomatic)

	# check the absence of _R1_P.fastq.gz and _R2_P.fastq.gz files of the step 3_trimming, then return warning or successfull message with the function absentefile_warning_success of the module genomics.py
	genomic.absentefile_warning_success(expectedfile = r + '/' + '3_trimming' + '/' + R1id + '_R1_P.fastq.gz')
	genomic.absentefile_warning_success(expectedfile = r + '/' + '3_trimming' + '/' + R1id + '_R2_P.fastq.gz')

	# prepare and run BWA
	bwaoutput = mappingoutput + R1id + '.sam'
	REname = os.path.basename(RE)
	cmdBWA = BW + ' mem' + ' -t ' + str(t) + ' -K 10000000 ' + referenceoutput + REname + ' ' + R1P + ' ' + R2P + ' > ' + bwaoutput
	print("#### Execute %s" %cmdBWA)
	os.system(cmdBWA)

	# check the absence of .sam file of the step 4_mapping, then return warning or successfull message with the function absentefile_warning_success of the module genomics.py
	genomic.absentefile_warning_success(expectedfile = r + '/' + '4_mapping' + '/' + R1id + '.sam')

	# prepare and run SortSam from Picard
	sortedbamoutput = mappingoutput + R1id + '.sorted.bam'
	cmdSortSam = PI + ' SortSam ' + '-I ' + bwaoutput + ' -O ' + sortedbamoutput + ' --SORT_ORDER coordinate' + ' --VALIDATION_STRINGENCY LENIENT' + ' --CREATE_INDEX true'
	print("#### Execute %s" %cmdSortSam)
	os.system(cmdSortSam)

	# check the absence of the .sorted.bam file of the step 4_mapping, then return warning or successfull message with the function absentefile_warning_success of the module genomics.py
	genomic.absentefile_warning_success(expectedfile = r + '/' + '4_mapping' + '/' + R1id + '.sorted.bam')

	# prepare and run "add or replace read group" (ARGG) from Picard
	arrgbamoutput = mappingoutput + R1id + '.ARRG.bam'
	cmdARRG = PI + ' AddOrReplaceReadGroups ' + '-I ' + sortedbamoutput + ' -O ' + arrgbamoutput + ' --SORT_ORDER coordinate' + ' --VALIDATION_STRINGENCY LENIENT' + ' --CREATE_INDEX true' + ' --RGID ' + R1id + ' --RGLB ' + R1id + ' --RGPU ' + R1id + ' --RGPL illumina' + ' --RGSM ' + R1id
	print("#### Execute %s" %cmdARRG)
	os.system(cmdARRG)

	# check the absence of the .ARGG.bam file of the step 4_mapping, then return warning or successfull message with the function absentefile_warning_success of the module genomics.py
	genomic.absentefile_warning_success(expectedfile = r + '/' + '4_mapping' + '/' + R1id + '.ARRG.bam')

	# prepare and run MarkDuplicates (MD) from Picard
	mdbamoutput = mappingoutput + R1id + '.dedupped.bam'
	metricsmdbamoutput = mappingoutput + R1id + '.dedupped.metrics'
	cmdMD = PI + ' MarkDuplicates ' + '-I ' + arrgbamoutput + ' -O ' + mdbamoutput + ' -M ' + metricsmdbamoutput + ' --CREATE_INDEX true' + ' --VALIDATION_STRINGENCY LENIENT' + ' --REMOVE_SEQUENCING_DUPLICATES true' + ' --REMOVE_DUPLICATES true'
	print("#### Execute %s" %cmdMD)
	os.system(cmdMD)

	# check the absence of the .dedupped.bam file of the step 4_mapping, then return warning or successfull message with the function absentefile_warning_success of the module genomics.py
	genomic.absentefile_warning_success(expectedfile = r + '/' + '4_mapping' + '/' + R1id + '.dedupped.bam')

	# prepare and run HaplotypeCaller (HA) from GATK4
	gcvfoutput = callingoutput + R1id + '.g.vcf.gz'
	cmdHC = GA + ' --java-options "-Xms50g -Xmx50g" HaplotypeCaller ' + '-R ' + referenceoutput + REname + ' -I ' + mdbamoutput + ' -O ' + gcvfoutput + ' -ERC GVCF -stand-call-conf 30 -ploidy 1 -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation'
	print("#### Execute %s" %cmdHC)
	os.system(cmdHC)

	# check the absence of the .g.vcf.gz file of the step 5_calling, then return warning or successfull message with the function absentefile_warning_success of the module genomics.py
	genomic.absentefile_warning_success(expectedfile = r + '/' + '5_calling' + '/' + R1id + '.g.vcf.gz')

	#calculate average depth coverage from the .dedupped.bam file
	depthout = qualityoutput + R1id + '.bam.depth'
	cmddepth = SA + ' depth ' + mdbamoutput + ''' |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average Depth Coverage (X) =",sum/NR; print "Standard Deviation Depth Coverage (X) =",sqrt(sumsq/NR - (sum/NR)**2)}' > ''' + depthout
	print("#### Execute %s" %cmddepth)
	os.system(cmddepth)

	#calculate breadth coverage from the .dedupped.bam file
	breadthout = qualityoutput + R1id + '.bam.metrics'
	cmdbreadth = SA + ' coverage ' + mdbamoutput + ' > ' + breadthout
	print("#### Execute %s" %cmdbreadth)
	os.system(cmdbreadth)

	#modify breadth coverage from the .bam.metrics file
	modifbreadthout = qualityoutput + R1id + '.bam.breadth'
	cmdmodifbreadth = 'cat ' + breadthout + ''' | awk ' { print $6 } ' | tr -d "\n" | sed 's@coverage@Breadth Coverage (%) = @' > ''' + modifbreadthout
	print("#### Execute %s" %cmdmodifbreadth)
	os.system(cmdmodifbreadth)

	#remove the .bam.metrics file
	cmdrmbreadthout = 'rm ' + breadthout
	print("#### Execute %s" %cmdrmbreadthout)
	os.system(cmdrmbreadthout)

	# congtratulate users with the function congratulation of the module genomic.py
	genomic.congratulation()

# driver code: if the code above is a scrypt, call  main() function, rather than to considere it as a module
if __name__ == "__main__":
	# calling main() function
	main()
