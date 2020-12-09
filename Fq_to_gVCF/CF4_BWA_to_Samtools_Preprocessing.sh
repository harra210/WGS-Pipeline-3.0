#!/bin/bash
##### USER INTERACTIVE AREA #####
# This prompts the user to verify that files that have been concatenated properly have had their names set properly
echo "Fastq naming convention:"
echo "Lab generated sequence: BreedName#_SPOTID"
echo "SRA or other lab generated sequence: BreedName#"
#
# Give the script a slight pause for the user to read the statements
#sleep 1
#
read -sp "`echo -e 'Are the fastq names set properly? If not press Ctrl+C. If so press any key to continue.\n\b'`" -n1 key
#
#This prompts the user to input whether the samples that are being put into the pipeline were ALL generated at NISC or known PCR-free library preparation.
#
#read -p "Are the samples being processed ALL generated at NISC or known to have been prepared using PCR-free libraries? (enter yes or no) " indel_modelQ
#while true; do
#	case "$indel_modelQ" in
#        	[YyEeSs]* ) indel_model="NONE"; break;;
#        	[NnOo]* ) indel_model="CONSERVATIVE"; break;;
#        	*) echo "Enter yes or no" ;;
#	esac
#done
#
# Prompts the user the input the parent directory where the concatenated fastq files are located with a 30 second timeout.
echo "What parent directory are your fastq files that you want to align?"
read -e -t 30 FQ_DIR
echo "What do you want to name your base swarm?"
read -e -t 30 SWARM_NAME
#
##### END OF USER INTERACTIVE SECTION ######
#
# Variable defining sections
homedir=$(pwd)
> bwa_to_picard.swarm
cd ../
basedir=$(pwd)
cd tmp/
tmpdir=$(pwd)
#cd RG
#RGdir=$(pwd)
CanFam31="/data/Ostrander/Resources/cf31PMc.fa"
CanFam4="/data/Ostrander/Resources/CanFam4_Tasha/Dog10K_Boxer_Tasha_1.0.KP081776.1.fa"
#
# In the tmp folder, blanking the used tmp files
> raw_files.tmp
> directories.tmp
> LB.tmp
> samplenames.tmp
> single_dir.tmp
> IDs.tmp
> filenames.tmp
> PU.tmp
# Change to the fastq directory to populate RG header files
#FQ DEBUG
#FQ_DIR="/data/harrisac2/Heidi_TCC/Fastq/Emergency/Test"
#
cd $FQ_DIR
# First find section takes care of the concatenated files. The directory search is shared across the script.
#find . -maxdepth 2 -name "*_1.fastq.gz" -printf '%f\n' | sed 's/_1.fastq.gz//' &> samplename.txt
find $PWD -maxdepth 2 -name "*L00?_R1_001.fastq.gz" -printf '%h\n' &> "$tmpdir"/directories.tmp
find $PWD -maxdepth 2 -name "*L001_R1_001.fastq.gz" -printf '%h\n' &> "$tmpdir"/single_dir.tmp
# Second find section will search the raw files and set them into their respective temp files and then parse the raw_files temp file to create the library file.
#
#
#This section creates the files necessary to input the bam headers later on.
#
cd $tmpdir
#
IFS=,$'\n' read -d '' -r -a singledir < single_dir.tmp
IFS=,$'\n' read -d '' -r -a multidir < directories.tmp
#
#
#DEBUG AREA
#echo ${directory[@]}
#sleep 30
for ((i = 0; i < ${#singledir[@]}; i++))
do
        cd ${singledir[$i]}; find . -maxdepth 2 -name "*L00?_R1_001*" -printf '%f\n' | sort -V >> "$tmpdir"/raw_files.tmp
done
#
#for ((i = 0; i < ${#singledir[@]}; i++))
#do
#	cd ${singledir[$i]}; > samplename.txt; find . -maxdepth 1 -name "*_1.fastq.gz" -printf '%f\n' | sed 's/_1.fastq.gz//' &> samplename.txt
#done
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
	cd ${singledir[$i]}; > raw_header.txt; find . -maxdepth 2 -name "*L00?_R1_001*" -printf '%f\n' | sed 's/_R1_001.fastq.gz//' | sort -V &> filenames.txt
done
#
awk -F '_' '{print $2}' "$tmpdir"/raw_files.tmp > "$tmpdir"/LB.tmp
#
cd $tmpdir
IFS=,$'\n' read -d '' -r -a LB < LB.tmp
IFS=,$'\n' read -d '' -r -a rawname < raw_files.tmp
#IFS=,$'\n' read -d '' -r -a samplename < samplenames.tmp
#
#for ((i = 0; i < ${#multidir[@]}; i++))
#do
#	echo ""${multidir[$i]}":"${samplename[$i]}""
#done
#sleep 30
#
#DEBUG
#echo ${LB[@]}
#echo ${rawname[@]}
#echo ${samplename[@]}
#sleep 30
# Iterate through each directory with fastqs and head the first line of the raw R1 fastq and pipe that into file named raw_header.txt that is placed in the samples folder
#
for ((i = 0; i < ${#multidir[@]}; i++))
do
	cd "${multidir[$i]}"; zcat "${rawname[$i]}" | head -n 1 >> "${multidir[$i]}"/raw_header.txt
#	cd "${directory[$i]}"; gunzip -c "${rawname[$i]}" | awk ' NR==1 {print; exit}' >> "${directory[$i]}"/raw_header.txt
done
#
## DEBUG SECTION
#for ((i = 0; i < ${#multidir[@]}; i++))
#do
#	cd "${multidir[$i]}"; gunzip -c "${rawname[$i]}" | awk ' NR==1 {print; exit}'
#done
#sleep 30
#
# Iterates through each directory again, and then out of that print out just the Flowcell and Sample Tag information and label that as the trimmed_header
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
        cd ${singledir[$i]}; awk -F ':' '{print $3,$4,$10}' raw_header.txt > trimmed_header.txt
done
#
# Iterate through the directories, and then in each directory create a file labeled LB_header.txt that because the arrays are linked due to the temporary file creation is correct and can be verified manually looking at the original filenames
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
        cd ${singledir[$i]}; printf "${LB[$i]}" > LB_header.txt
done
#
# Iterate through the directories and pasting the two files trimmed_header and LB_header into the pasted_header.txt which will have all of the Read Group information that
for ((i = 0; i < ${#singledir[@]}; i++))
do
        cd ${singledir[$i]}; awk '{print $1"."$2"."$3}' trimmed_header.txt > final_PU.txt
done
#
#IFS=,$'\n' read -d '' -r -a PU < final_PU.txt
#IFS=,$'\n' read -d '' -r -a filename < filenames.txt
#
for ((i = 0; i < ${#multidir[@]}; i++))
do
	while read -r SM
	do
	cd "${multidir[$i]}"; echo "$SM" >> "$tmpdir"/IDs.tmp
	done < "${multidir[$i]}"/samplename.txt
done
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
	cd "${singledir[$i]}"; cat final_PU.txt >> "$tmpdir"/PU.tmp
done
#
for ((i = 0; i < ${#singledir[@]}; i++))
do
	cd "${singledir[$i]}"; cat filenames.txt >> "$tmpdir"/filenames.tmp
done
#sleep 30
IFS=,$'\n' read -d '' -r -a ID < "$tmpdir"/IDs.tmp
IFS=,$'\n' read -d '' -r -a PU < "$tmpdir"/PU.tmp
IFS=,$'\n' read -d '' -r -a samplename < "$tmpdir"/filenames.tmp
#bamname=("${samplename[@]}")
#declare -a bamname
#
#for ((i = 0; i < "${#bamname[@]}"; i++))
#do
#	echo "${bamname[$i]}"_"${ID[$i]}"
#done
#echo ${ID[@]}
#echo ${samplename[@]}
#sleep 15
#unset bamname
#for (( i=0; i<${#samplename[*]}; ++i));
#do
#	bamname+=( "${ID[$i]}""_""${samplename[$i]}" )
#done
#echo ${bamname[@]}
#echo ${!bamname[@]}
#echo ${!samplename[@]}
#echo ${!ID[@]}
#sleep 30
#This finishes the creating of proper arrays.
#
#This section begins creating the swarmfiles for processing.
#
cd $homedir
#
for ((i = 0; i < ${#multidir[@]}; i++))
do
	echo "cd ${multidir[$i]}; bwa mem -M -t \$SLURM_CPUS_PER_TASK -R '@RG\\tID:"${ID[$i]}"\\tSO:coordinate\\tLB:"${LB[$i]}"\\tPL:ILLUMINA\\tSM:"${ID[$i]}"\\tPU:"${PU[$i]}"' "$CanFam4" "${samplename[$i]}"_R1_001.fastq.gz "${samplename[$i]}"_R2_001.fastq.gz | samtools view -h | samtools sort -@ \$SLURM_CPUS_PER_TASK -T /lscratch/\$SLURM_JOB_ID/${samplename[$i]} -o "${multidir[$i]}"/sort_"${ID[$i]}"_"${samplename[$i]}".bam && samtools flagstat "${multidir[$i]}"/sort_"${ID[$i]}"_"${samplename[$i]}".bam > sort_"${ID[$i]}"_"${samplename[$i]}".flagstat" >> "$homedir"/bwa_to_picard.swarm
done
#
cd $homedir
#
head bwa_to_picard.swarm
read -p "Does the formatting of the swarmfile appear correct? (yes or no) " promptA
while true; do
        case "$promptA" in
                [YyEeSs]* ) break ;;
                [NnOo]* ) echo "Verify that the inputs are correct and try again"; exit;;
                *) echo "Enter yes or no" ;;
        esac
done
#
echo "Swarm Job ID: "
#
jobid1=$(swarm -f bwa_to_picard.swarm -g 36 -t 20 --gres=lscratch:350 --time 2-0 --module bwa,samtools,picard,GATK/4.1.9.0 --logdir ~/job_outputs/bwa_to_picard/"$SWARM_NAME"_FQ --sbatch "--mail-type=ALL,TIME_LIMIT_90 --job-name "$SWARM_NAME"_FQ")
echo $jobid1
#
