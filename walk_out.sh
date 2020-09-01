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
#Variable setting
pwd=$(pwd)
REF="/data/Ostrander/Resources/"
IDX="/data/Ostrander/Resources/"
FQDIR="N/A"
echo "Where do you want to output your files?"
read -ep "Output directory: " MAPDIR
cd $MAPDIR
mkdir -p Logs
cd Logs
LOGDIR=$(pwd)
cd $MAPDIR
mkdir -p Metrics
cd Metrics
MEDIR=$(pwd)
cd $MAPDIR
mkdir -p Qual
cd Qual
QUALDIR=$(pwd)
cd $pwd
#
> out.file
bash walk_test.sh | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]//g" > out.file
sed -i '$ d' out.file
#
cd tmp/
tmpdir=$(pwd)
> SM.var
> IDfl.var
> IDlb.var
> IDln.var
> Rd.var
cd filldown
fill=$(pwd)
echo "$REF" &> 1_REF.var
echo "$IDX" &> 2_IDX.var
echo "$FQDIR" &> 3_FQDIR.var
echo "$MAPDIR" &> 4_MAPDIR.var
echo "$LOGDIR" &> 5_LOGDIR.var
echo "$MEDIR" &> 6_MEDIR.var
echo "$QUALDIR" &> 7_QUALDIR.var
#
find $PWD -name "*.var" | sort -V &> "$tmpdir"/filldown.txt
#
cd ../../
homedir=$(pwd)
#
awk -F '\t' '{print $1}' out.file > "$tmpdir"/SM.var
awk -F '_' '{print $1}' out.file > "$tmpdir"/IDfl.var
awk -F '_' '{print $2}' out.file > "$tmpdir"/IDlb.var
awk -F '_' '{print $4}' out.file | sed 's\L00\\' > "$tmpdir"/IDln.var
awk -F '_' '{print $5}' out.file | sed 's\R\\' > "$tmpdir"/Rd.var
awk -F '\t' '{print $3}' out.file > "$tmpdir"/D1.var
#

awk -F '\t' '{print $2}' out.file > "$tmpdir"/D2.var
#
cd $tmpdir
paste IDfl.var IDln.var IDlb.var SM.var Rd.var D1.var D2.var > "$homedir"/paste_test.txt
#
cd $homedir
#cut -d $'\t' -f 1 paste_test.txt > paste_test_cut.txt
awk -F '\t' '{print $2,$3,$4,$5,$6,$7,$8}' paste_test.txt > paste_test_cut.txt
#
cp paste_test_cut.txt paste_test_add.txt
#
cd $tmpdir
#
#This part will head the filldown portion and then shuttle the elements of the filldown files into an array which can then be used to paste the last half
files=( $(while IFS= read -r line; do head -n 1 "$line"; done < filldown.txt))
declare -a files
#
#echo ${files[@]}
#
cd $homedir
for f in ${files[@]}
do
        printf "%s " "$f" >> filldown_top.txt
done
#
cpRow=$(cat paste_test_add.txt | wc -l)
#echo $cpRow
#more files_test.txt
#printf "%s\n" ${files[@]} >> paste_test_add.txt
> final.txt
for x in $cpRow
do
        awk '{for (i=1; i<=n; i++) print}' n="$x" files_test.txt &> files_test_awk.txt
done
#more files_test_awk.txt
#sleep 10
#
paste paste_test_add.txt files_test_awk.txt > final.txt
sed -i $'s/\t/ /g' final.txt
#sed -i $'s/ /,/g' final.txt
#more final.txt
sed -i '1i IDfl IDln IDlb SM R1 R2 D1 D2 REF IDX FQDIR MAPDIR LOGDIR MEDIR QUALDIR' final.txt
