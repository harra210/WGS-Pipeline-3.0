#!/bin/bash
# Setting initial variables and blanking temporary files
#
homedir=$(pwd)
> files.tmp
> directories.tmp
> sample_swarmfile.txt
#
# Current hard coded location for Fastqs, will change to softcode in final version
#
cd /data/harrisac2/Heidi_TCC/Fastq/Emergency
#
# Change directory and set the parent directory for fastq location
#
fastqdir=$(pwd)
#
# Use Linux find command for the first lane, read 1 of fastqs and pipe that to temporary files
find . -maxdepth 2 -name "*L001_R1_001*" -printf '%f\n' &> $homedir/files.tmp
find $PWD -maxdepth 2 -name "*L001_R1_001*" -printf '%h\n' &> $homedir/directories.tmp
#
cd $homedir
#
# Take column two from the results of the file search piping to a new tmp file with Library information
#
awk -F '_' '{print $2}' files.tmp > LB.tmp
#
# Setting contents of the temporary files into respective arrays
IFS=,$'\n' read -d '' -r -a filename < files.tmp
IFS=,$'\n' read -d '' -r -a directory < directories.tmp
IFS=,$'\n' read -d '' -r -a LB < LB.tmp
#
# Iterate through each directory with fastqs and head the first line of the raw R1 fastq and pipe that into file named raw_header.txt that is placed in the samples folder
#
for ((i = 0; i < ${#directory[@]}; i++))
do
	cd ${directory[$i]}; zcat ${filename[$i]} | head -n 1 &> raw_header.txt
done
#
# Iterates through each directory again, and then out of that print out just the Flowcell and Sample Tag information and label that as the trimmed_header
#
for ((i = 0; i < ${#directory[@]}; i++))
do
	cd ${directory[$i]}; awk -F ':' '{print $3,$10}' raw_header.txt > trimmed_header.txt
done
#
# Iterate through the directories, and then in each directory create a file labeled LB_header.txt that because the arrays are linked due to the temporary file creation is correct and can be verified manually looking at the original filenames
#
for ((i = 0; i < ${#directory[@]}; i++))
do
	cd ${directory[$i]}; printf "${LB[$i]}" > LB_header.txt
done
#
# Iterate through the directories and pasting the two files trimmed_header and LB_header into the pasted_header.txt which will have all of the Read Group information that 
for ((i = 0; i < ${#directory[@]}; i++))
do
	cd ${directory[$i]}; paste -d ' ' trimmed_header.txt LB_header.txt > pasted_header.txt
done
#
# Final iteration for each directory iterated change into it and then read the pasted_header.txt in that directory and set temporary variables as FC (Flowcell), TAG (Sample Tag), and LB (Library) and then echo a string as a sample swarmfile containing the Read Group information for each sample.
#
for ((i = 0; i < ${#directory[@]}; i++))
do
	cd ${directory[$i]};
	while IFS=" " read -r FC TAG LB
	do
		echo "cd "${directory[$i]}"; RGID="$FC" RGPU="$FC"."$TAG" RGLB="$LB"" >> "$homedir"/sample_swarmfile.txt
	done < "${directory[$i]}"/pasted_header.txt
done
