#!/bin/bash
##### USER INTERACTIVE AREA #####
# This prompts the user to verify that files that have been concatenated properly have had their names set properly
echo "Fastq naming convention:"
echo "Lab generated sequence: BreedName#_SPOTID"
echo "SRA or other lab generated sequence: BreedName#"
#
# Give the script a slight pause for the user to read the statements
sleep 1
#
read -sp "`echo -e 'Are the fastq names set properly? If not press Ctrl+C. If so press any key to continue.\n\b'`" -n1 key
#
#This prompts the user to input whether the samples that are being put into the pipeline were ALL generated at NISC or known PCR-free library preparation.
#
read -p "Are the samples being processed ALL generated at NISC or known to have been prepared using PCR-free libraries? (enter yes or no) " indel_modelQ
while true; do
	case "$indel_modelQ" in
        	[YyEeSs]* ) indel_model="NONE"; break;;
        	[NnOo]* ) indel_model="CONSERVATIVE"; break;;
        	*) echo "Enter yes or no" ;;
	esac
done
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
#
# In the tmp folder, blanking the used tmp files
> raw_files.tmp
> concat_fq_files.tmp
> directories.tmp
> LB.tmp
# Change to the fastq directory to populate RG header files
#
cd $FQ_DIR
# First find section takes care of the concatenated files. The directory search is shared across the script.
find . -maxdepth 2 -name "*_1.fastq.gz" -printf '%f\n' | sed 's/_1.fastq.gz//' &> "$tmpdir"/concat_fq_files.tmp
find $PWD -maxdepth 2 -name "*_1.fastq.gz" -printf '%h\n' &> "$tmpdir"/directories.tmp
# Second find section will search the raw files and set them into their respective temp files and then parse the raw_files temp file to create the library file.
#
#
#This section creates the files necessary to input the bam headers later on.
#
cd $tmpdir
#
IFS=,$'\n' read -d '' -r -a directory < directories.tmp
#
for ((i = 0; i < ${#directory[@]}; i++))
do
        cd ${directory[$i]}; find . -maxdepth 2 -name "*L001_R1_001*" -printf '%f\n' -quit >> "$tmpdir"/raw_files.tmp
done
#
awk -F '_' '{print $2}' "$tmpdir"/raw_files.tmp > "$tmpdir"/LB.tmp
#
cd $tmpdir
IFS=,$'\n' read -d '' -r -a LB < LB.tmp
IFS=,$'\n' read -d '' -r -a rawname < raw_files.tmp
IFS=,$'\n' read -d '' -r -a filename < concat_fq_files.tmp
#
# Iterate through each directory with fastqs and head the first line of the raw R1 fastq and pipe that into file named raw_header.txt that is placed in the samples folder
#
for ((i = 0; i < ${#directory[@]}; i++))
do
        cd "${directory[$i]}"; zcat "${rawname[$i]}" | head -n 1 &> "${directory[$i]}"/raw_header.txt
done
#
## DEBUG SECTION
#for ((i = 0; i < ${#directory[@]}; i++))
#do
#	echo "${directory[$i]}; zcat "${rawname[$i]}" | head -n 1"
#done
#sleep 30
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
        cd ${directory[$i]}; paste -d ' ' trimmed_header.txt LB_header.txt > final_header_variables.txt
done
#
#This finishes the creating of proper arrays.
#
#This section begins creating the swarmfiles for processing.
#
cd $homedir
#
for ((i = 0; i < ${#directory[@]}; i++))
do
	cd ${directory[$i]};
	while IFS=" " read -r FC TAG LB
	do
		echo "cd ${directory[$i]}; bwa mem -M -t \$SLURM_CPUS_PER_TASK "$CanFam31" "${filename[$i]}"_1.fastq.gz "${filename[$i]}"_2.fastq.gz | samtools view -h | samtools sort -@ \$SLURM_CPUS_PER_TASK -T /lscratch/\$SLURM_JOB_ID/${filename[$i]} -o /lscratch/\$SLURM_JOB_ID/sort_"${filename[$i]}".bam && samtools flagstat /lscratch/\$SLURM_JOB_ID/sort_"${filename[$i]}".bam > "${filename[$i]}".flagstat && java -Xmx4g -jar \$PICARDJARPATH/picard.jar CollectMultipleMetrics I=/lscratch/\$SLURM_JOB_ID/sort_"${filename[$i]}".bam O="${filename[$i]}".AlignmentMetrics R="$CanFam31" && java -Xmx16g -jar \$PICARDJARPATH/picard.jar AddOrReplaceReadGroups I=/lscratch/\$SLURM_JOB_ID/sort_"${filename[$i]}".bam O=/lscratch/\$SLURM_JOB_ID/RG_"${filename[$i]}".bam SO=coordinate RG="${filename[$i]}" RGLB="${LB[$i]}" RGPL=ILLUMINA RGSM="${filename[$i]}" RGPU="$FC"."$TAG" && rm /lscratch/\$SLURM_JOB_ID/sort_"${filename[$i]}".bam && java -Xmx16g -jar \$PICARDJARPATH/picard.jar MarkDuplicates I=/lscratch/\$SLURM_JOB_ID/RG_"${filename[$i]}".bam O="${directory[$i]}"/dedup_"${filename[$i]}".bam M="${directory[$i]}"/"${filename[$i]}"_metrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true TMP_DIR=/lscratch/\$SLURM_JOB_ID && samtools index dedup_"${filename[$i]}".bam" >> "$homedir"/bwa_to_picard.swarm
	done < "${directory[$i]}"/final_header_variables.txt
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
#jobid1=$(swarm -f bwa_to_picard.swarm -g 36 -t 20 --gres=lscratch:350 --time 2-0 --module bwa, samtools, picard --logdir ~/job_outputs/bwa_to_picard/"$SWARM_NAME"_FQ --sbatch "--mail-type=ALL, TIME_LIMIT_90 --job-name "SWARM_NAME"_FQ")
#echo $jobid1
#
#cd gatk/BaseRecalibrator
#gatkBRdir=$(pwd)
> gatk_BRPR.swarm
#
for ((i = 0; i < ${#directory[@]}; i++))
do
	echo "cd ${directory[$i]}/; gatk --java-options \"-Xmx16g -XX:ParallelGCThreads=4\" BaseRecalibrator -bqsr-baq-gap-open-penalty 30.0 -R "$CanFam31" --tmp-dir /lscratch/\$SLURM_JOB_ID -I dedup_"${filename[$i]}".bam --known-sites /data/Ostrander/Resources/CFA31_151.dbSNP_numb_order.vcf -O "${filename[$i]}"_recal.table; gatk --java-options \"-Xmx16g -XX:ParallelGCThreads=4\" ApplyBQSR -R "$CanFam31" --tmp-dir /lscratch/\$SLURM_JOB_ID -I dedup_"${filename[$i]}".bam -bqsr "${filename[$i]}"_recal.table -O "${filename[$i]}"_BQSR.bam -OBM && rm dedup_"${filename[$i]}".bam && rm dedup_"${filename[$i]}".bam.bai; samtools depth "${filename[$i]}"_BQSR.bam | awk '{sum+=\$3} END {print sum/NR}' > "${filename[$i]}".coverageALL; samtools depth -r chrX "${filename[$i]}"_BQSR.bam | awk '{sum+=\$3| END {print sum/NR}' > "${filename[$i]}".coveragechrX" >> gatk_BRPR.swarm
done
head gatk_BRPR.swarm
read -sp "`echo -e 'Verify syntax of the Basecalibrator swarm file \n\b'`" -n1 key
#echo "BaseRecalibrator Swarm JobID:"
#jobid2=$(swarm -f gatk_BRPR.swarm -g 18 -t 6 --time 4-0 --gres=lscratch:50 --module samtools,GATK/4.1.8.0 --logdir ~/job_outputs/gatk/BaseRecalibrator/"$SWARM_NAME"_BQSR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_BQSR --dependency=afterok:$jobid1")
#echo $jobid2
#
# HaplotypeCaller Section
> gatk_HCaller.swarm
for ((i = 0; i < ${#directory[@]}; i++))
do
	echo "cd ${directory[$i]}; gatk --java-options \"-Xmx12g\" HaplotypeCaller -R "$CanFam31" -I "${filename[$i]}"_BQSR.bam -O "${filename[$i]}"_g.vcf.gz --output-mode EMIT_ALL_ACTIVE_SITES -ERC GVCF --pcr-indel-model $indel_model --smith-waterman FASTEST_AVAILABLE --tmp-dir /lscratch/\$SLURM_JOB_ID -OVM" >> gatk_HCaller.swarm
done
head gatk_HCaller.swarm
