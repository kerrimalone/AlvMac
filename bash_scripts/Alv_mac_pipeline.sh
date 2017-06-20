##########################################################################
# RNA-seq Time Course: M. tuberculosis V M. bovis infected bovine		 #
# alveolar macrophages 												 #
# paired-end reads.                                                      #
#       --- Linux bioinformatics workflow for known sense genes ---      #
##########################################################################
# Based on the pipeline created by Nalpas, N.C. (2014) 
# DOI badge: http://dx.doi.org/10.5281/zenodo.12474
# Author of current version (1.0.0): Malone, K.M


# All files were downloaded from MSU in 2013. At the time, they were
# also md5sum checked and renamed.
#Copied the raw sequencing files from 3 sequencing pools from 
# /home/nnalpas/storage/ALV_MAC_RNAseq/fastq_sequence/raw_reads
#Copy the 3 pool adapter sequence files from dmagee
cp /home/dmagee/scratch/ALV_MAC_RNAseq/fastq_sequence/pool_*.txt ~/storage/RNAseq/Alv_mac_data/raw_reads/.


###########################################
# FastQC quality check of raw FASTQ files #
###########################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/poolA
cd !$ #switch into the newly created directory

# Create and enter the quality check output directory:
mkdir -p $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/poolB
cd !$ #switch into the newly created directory

# Create and enter the quality check output directory:
mkdir -p $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/poolC
cd !$ #switch into the newly created directory

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/poolA \
--noextract --nogroup -t 2 \
$HOME/storage/RNAseq/Alv_mac_data/raw_reads/DM-pool-A/Raw_FCC1VDGACXX-CHKPEI12120043_L1_1.fq.gz

#Check output and proceed if acceptable.

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find $HOME/storage/RNAseq/Alv_mac_data/raw_reads/DM-pool-A/ \
-name *fq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/poolA $file" \
>> fastqc_poolA.sh; done;

for file in `find $HOME/storage/RNAseq/Alv_mac_data/raw_reads/DM-pool-B/ \
-name *fq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/poolB $file" \
>> fastqc_poolB.sh; done;

for file in `find $HOME/storage/RNAseq/Alv_mac_data/raw_reads/DM-pool-C/ \
-name *fq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/poolC $file" \
>> fastqc_poolC.sh; done;

# Split and run all scripts for poolA:
split -d -l 70 fastqc_poolA.sh fastqc_poolA.sh.
for script in `ls fastqc_poolA.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Split and run all scripts for poolB:
split -d -l 70 fastqc_poolB.sh fastqc_poolB.sh.
for script in `ls fastqc_poolB.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Split and run all scripts for poolC:
split -d -l 70 fastqc_poolC.sh fastqc_poolC.sh.
for script in `ls fastqc_poolC.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all the files were processed:
for file in `ls fastqc_poolA.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc_poolA.txt
done

# Check if all the files were processed:
for file in `ls fastqc_poolB.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc_poolB.txt
done

# Check if all the files were processed:
for file in `ls fastqc_poolC.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc_poolC.txt
done

# Deleted all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/tmp

for file in `find ~/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/poolA/*_fastqc.zip`; do unzip \
$file -d $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/tmp; \
done;

for file in \
`find $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; done

for file in \
`find $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Remove temporary folder and its files:
rm -rf $HOME/scratch/RNAseq/Alv_mac/quality_check/pre-filtering/tmp


##################################################################
# Deconvolution of sequencing reads using FASTX Barcode splitter #
##################################################################
#Example (Assuming 's_2_100.txt' is a FASTQ file, 'mybarcodes.txt' is the barcodes file):            
#$ cat s_2_100.txt | /usr/local/bin/fastx_barcode_splitter.pl \
#--bcfile mybarcodes.txt \
#--bol \
#--mismatches 2 \
#--prefix /tmp/bla_ --suffix ".txt"

# index files copied from dmagee above are in wrong format, columns need to be swapped.
# awk '{ print $2 "\t " $1}' pool_A_indices.txt > poolA_indices.txt

# #Create and enter the quality check output directory:
# mkdir -p $HOME/scratch/RNAseq/Alv_mac/deconv/poolA
# cd !$ #switch into the newly created directory

# #test the code on one file to make sure it is working
# zcat Raw_FCC1VDGACXX-CHKPEI12120043_L1_1.fq.gz | /usr/local/bin/fastx_barcode_splitter.pl \
# --bcfile $HOME/storage/RNAseq/Alv_mac_data/raw_reads/poolA_indices.txt \
# --bol \
# --mismatches 2 \
# --prefix $HOME/scratch/RNAseq/Alv_mac/deconv/poolA \
# --suffix ".fq.gz"


##############################################################################
# Alignment of FASTQ files against the Bos taurus reference genome with STAR #
##############################################################################

# Required software is STAR 2.5.1b, consult manual/tutorial for details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Download Bos taurus reference genome, version UMD3.1.1 from NCBI:
mkdir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/source_file
cd /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/source_file
nohup wget -o logfile -r -nd \
"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz" \
&
gunzip GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz

# Download annotation file for UMD3.1.1 NCBI Bos taurus Annotation Release 105:
mkdir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file
cd /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz
gunzip GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz

# Generate genome indexes files using annotations:
mkdir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index
cd /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index

nohup STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index \
--genomeFastaFiles \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/source_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna \
--sjdbGTFfile /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 \
--outFileNamePrefix \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index/Btau-UMD3.1.1 &

# Create and enter alignment working directory:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment
cd $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment


#STAR command for one pair of files, check it works by looking at % mapping in log.out
nohup STAR --runMode alignReads --runThreadN 10 --genomeLoad LoadAndRemove \
--genomeDir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index/ \
--readFilesIn \
$HOME/storage/RNAseq/Alv_mac_data/deconv_raw_reads/N1178_CN_24H_pe1.fastq.gz,\
$HOME/storage/RNAseq/Alv_mac_data/deconv_raw_reads/N1178_CN_24H_pe2.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./N1178_CN_24H_ \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx &


#run this block of code below to make sure you are pulling out the right files and creating
#correct folder name. names will be printed to screen.
for file in `find $HOME/storage/RNAseq/Alv_mac_data/deconv_raw_reads \
-name *_pe1.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/(_pe.)/_pe2/'`; \
foldername=`basename $file | perl -p -e 's/\_pe.\.fastq\.gz/_alignment/'`; \
echo $file
echo $file2;
echo $foldername; 
done; 


for file in `find $HOME/storage/RNAseq/Alv_mac_data/deconv_raw_reads \
-name *_pe1.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/(_pe.)/_pe2/'`; \
foldername=`basename $file | perl -p -e 's/\_pe.\.fastq\.gz/_alignment/'`; \
echo "mkdir $HOME/scratch/RNAseq/Alv_mac/STAR-2.5.1b_alignment/$foldername; \
cd $HOME/scratch/RNAseq/Alv_mac/STAR-2.5.1b_alignment/$foldername; \
STAR --runMode alignReads --runThreadN 10 --genomeLoad LoadAndRemove \
--genomeDir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index/ \
--readFilesIn \
$file,$file2 \
--readFilesCommand gunzip -c \
--outFilterMultimapNmax 10 --outFilterMismatchNmax 10 \
--outFileNamePrefix ./${foldername}_ --outSAMtype BAM Unsorted \
--outSAMattrIHstart 0 --outSAMattributes Standard --outReadsUnmapped Fastx" \
>> alignment.sh; \
done;


split -d -l 127 alignment.sh alignment.sh.
#STAR uses ~30G memory to load genome so only run one script at a time and run
#next script once memory drops below 30G and mapping begins
chmod 755 alignment.sh
nohup ./alignment.sh > alignment.sh.nohup &
#second run
chmod 755 alignment.sh.00
nohup ./alignment.sh.00 > alignment.sh.00.nohup &


# Check nohup.out file to see how many jobs finished successfully, will equal no of libraries(124):
grep -c 'Finished successfully' alignment.sh.00.nohup

# Merge all STAR log.final.out files into a single file using NNalpas custom script to parse STAR results:
#Reports Sample_id, No. of input reads, uniquely mapped reads etc. See NNalpas Github for script.
for file in `find $HOME/scratch/RNAseq/Alv_mac/STAR-2.5.1b_alignment \
-name *Log.final.out`; \
do perl /home/nnalpas/SVN/star_report_opener.pl -report $file; done;

#############################################
# FastQC quality check of aligned BAM files #
#############################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and go to working directory:
mkdir $HOME/scratch/RNAseq/Alv_mac/quality_check/post_alignment
cd $HOME/scratch/RNAseq/Alc_mac/quality_check/post_alignment

# Create a bash script to perform FastQC quality check on aligned SAM files:
for file in `find $HOME/scratch/RNAseq/Alv_mac/STAR-2.5.1b_alignment \
-name *.bam`; do echo "fastqc --noextract --nogroup -t 2 \
-o $HOME/scratch/RNAseq/Alv_mac/quality_check/post_alignment $file" >> \
fastqc_aligned.sh; done;

# Split and run all scripts 
split -d -l 35 fastqc_aligned.sh fastqc_aligned.sh.
for script in `ls fastqc_aligned.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Delete all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/scratch/RNAseq/Alv_mac/quality_check/post_alignment/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/RNAseq/Alv_mac/quality_check/post_alignment/tmp; \
done

for file in \
`find $HOME/scratch/RNAseq/Alv_mac/quality_check/post_alignment/tmp \
-name summary.txt`; do more $file >> reports_post-alignment.txt; done

for file in \
`find $HOME/scratch/RNAseq/Alv_mac/quality_check/post_alignment/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post_alignment.txt; \
done

# Check if all files were processed:
grep -c '##FastQC' basic_stats_post_alignment.txt
grep -c 'Basic Statistics' reports_post-alignment.txt
grep -c 'Analysis complete' fastqc_aligned.sh.00.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.01.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.02.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.03.nohup

# Remove temporary folder:
rm -r tmp/

###################################################################
# Summarisation of gene counts with featureCounts for sense genes #
###################################################################

# Required package is featureCounts, which is part of Subread 1.5.0-p1 software,
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf

# Create working directories:
cd $HOME/scratch/RNAseq/Alv_mac/
mkdir -p Count_summarisation/sense
cd $HOME/scratch/RNAseq/Alv_mac/Count_summarisation/sense

# Run featureCounts with one sample to check if it is working fine:
# -B, count those read pairs that have both ends aligned only
# -p, count read pairs rather than individual reads
# -C, Do not count read pairs that have their two ends mapping to diff chromos
# -R, Output detailed assignment
# -s 1, strandedness 
featureCounts -a \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
-B -p -C -R -s 1 -T 15 -t gene -g Dbxref -o ./counts.txt \
$HOME/scratch/RNAseq/Alv_mac/STAR-2.5.1b_alignment/N1178_CN_48H_alignment/N1178_CN_48H_alignment_Aligned.out.bam

# Create a bash script to run featureCounts on BAM file containing multihits and
# uniquely mapped reads using the reversely stranded parameter:
for file in `find $HOME/scratch/RNAseq/Alv_mac/STAR-2.5.1b_alignment \
-name *_Aligned.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Aligned.out.bam//'`; \
echo "mkdir $HOME/scratch/RNAseq/Alv_mac/Count_summarisation/sense/$sample; \
cd $HOME/scratch/RNAseq/Alv_mac/Count_summarisation/sense/$sample; \
featureCounts -a \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
-B -p -C -R -s 1 -T 10 -t gene -g Dbxref \
-o ${sample}_sense-counts.txt $file" >> sense_count.sh; done

# Split and run all scripts on Stampede:
split -d -l 70 sense_count.sh sense_count.sh.
for script in `ls sense_count.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
grep -c 'Read assignment finished.' sense_count.sh.00.nohup
grep -c 'Read assignment finished.' sense_count.sh.01.nohup


# Rename .featureCounts out files:
for folder in \
`ls $HOME/scratch/RNAseq/Alv_mac/Count_summarisation/sense/`; \
do file2=`echo $folder`; \
echo "cd $HOME/scratch/RNAseq/Alv_mac/Count_summarisation/sense/$file2; \
mv ./*_Aligned.out.bam.featureCounts ./${file2}_Aligned.sam.featureCounts" >> \
rename.sh; done;

chmod 755 rename.sh
nohup ./rename.sh > log.out &


# Create bash script to merge stats info from .featureCounts from all samples
# into a single file:
for file in `find $HOME/scratch/RNAseq/Alv_mac/Count_summarisation/sense/ \
-name *.featureCounts`; do echo echo \
"\`basename $file\` \`cut $file -f2 | sort | uniq -c | perl -p -e 's/\n/ /'\` >> \
annotation_summary_sense.txt" >> annotation_summary_sense.sh
done

# Split and run scripts on Stampede:
split -d -l 70 annotation_summary_sense.sh annotation_summary_sense.sh.
for script in `ls annotation_summary_sense.sh.*`
do
chmod 755 $script
nohup ./$script &
done

# Check that all files were processed:
grep -c '.featureCounts' annotation_summary_sense.txt

# Copy all *sense-counts.txt files to temporary folder:
mkdir $HOME/scratch/RNAseq/Alv_mac/Count_summarisation/sense/tmp

for file in `find $HOME/scratch/RNAseq/Alv_mac/Count_summarisation/sense/ \
-name *sense-counts.txt`; do cp $file \
-t $HOME/scratch/RNAseq/Alv_mac/Count_summarisation/sense/tmp; \
done

# Transfer all files for R analysis and then remove tmp folder:
rm -r tmp




























