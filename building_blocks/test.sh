#!/bin/bash
#
# This script is to generate a sample sheet that would be fitted to an individual sample which can then be fed into a pipeline and set up an array
# What this script will export should be a tab-delimited file with the following fields:
#
# IDfl: flowcell ID
# IDln: lane number
# SM: sample name
# R1: filename for read 1 (uncompressed fastq)
# R2: filename for read 2 (uncompressed fastq)
# D1: directory for read 1 / dir name for uBAM
# D2: directory for read 2 / file name for uBAM
# REF: ref genome base directory (note: there should be a target loci file included in this directory)
# IDX: Indexed reference genome directory
# FQDIR: directory for uncompressed fastq files
# MAPDIR: Output directory
# LOGDIR: log file directory
# MEDIR: Metrics file directory
# QUALDIR: directory for quality metrics
#
# Important note: The target loci file is a file that contains intervals that cover the entire genome, which will be important in how data is split up for a scatter gather approach in indel realignment.
#
#---------------------------------------------------
> files.txt
> SM.tmp
> SM.txt
START=$(date +%g-%m-%d_%H.%M.%S)
REF=/data/Ostrander/Resources/
IDX=/data/Ostrander/Resources/
LOGDIR=/data/harrisac2/Pipeline_CanFam4.0_Testing/
homedir=$(pwd)
#cd tmp
#tmpdir=$(pwd)
#cd $homedir
echo "What directory are your sample fastqs located in?"
read -ep "Fastq directory: " FILEDIR
#
cd $FILEDIR
# Set the the directories and take their names and then place them into an array for iterations
find . -type d -not -path "." -name "*" -printf '%f\n' &> $homedir/SM.tmp
cd $homedir
#sort -V SM.tmp > SM.txt
#$IFS=,$'\n' read -d '' -r -a list < SM.txt
#declare -a list
#Iterate across the arrays to set the variables
cd $FILEDIR
find . -name "*gz" -printf '%f\n' &> $homedir/files.txt
cd  $homedir
IFS=,$'\n' read -d '' -r -a files < files.txt
declare -a files
#
fl=$(awk -F '_' 'NR==1{print $1}' files.txt)
ln=$(awk -F '_' 'NR==1{print $4}' files.txt | sed 's\L00\\')
lb=$(awk -F '_' 'NR==1{print $2}' files.txt)
R=$(awk -F '_' 'NR==1{print $5}' files.txt | sed 's\R\\')
#echo $fl
#echo $ln
#echo $lb
#echo $R
#
#SAMPLE_SHEET="SM.txt"
#echo $SAMPLE_SHEET
#for ROW in $(seq 1 $(wc -l < $SAMPLE_SHEET))
#do
#	read IDfl IDln LB SM R1 R2 D1 D2 REF IDX FQDIR MAPDIR LOGDIR MEDIR QUALDIR <<< $(sed "${ROW}q;d" $SAMPLE_SHEET)
#done 
##
# Create a log directory so logs can go somewhere
#if [ -z ${START+x} ];
#	then
#		echo "START time is empty, check script cmd, exiting"
#		exit
#	else
#	fi
#

