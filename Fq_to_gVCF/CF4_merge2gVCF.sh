#!/bin/bash
echo "What parent directory are your preprocessed files located?"
read -e -t 30 FILE_DIR
echo "What do you want to name your base swarm?"
read -e -t 30 SWARM_NAME
#
homedir=$(pwd)
> CF4_mergeDedup.swarm
> CF4_BQSR.swarm
> CF4_HC.swarm
cd ../
basedir=$(pwd)
cd tmp/
tmpdir=$(pwd)
#
#Reference Genome Variables
CanFam31="/data/Ostrander/Resources/cf31PMc.fa"
CanFam4="/data/Ostrander/Resources/CanFam4_Tasha/Dog10K_Boxer_Tasha_1.0.KP081776.1.fa"
#
> merge_files.tmp
> merge_dir.tmp
> merge_IDs.tmp
> merge_lsfiles.tmp
> merge_samplename.tmp
> BaseRecalibrator_1pass.intervals
> BaseRecalibrator_2pass.intervals
#
#DEBUG $FILE_DIR
#FILE_DIR="/data/harrisac2/Heidi_TCC/Fastq/Emergency/Test"
#INT="/data/Ostrander/Alex/Intervals/CanFam4/interval_6.bed"
#INT_2="/data/Ostrander/Alex/Intervals/CanFam4/interval_35.bed"
#knownsite="/data/Ostrander/Resources/CanFam4_Tasha/Dog10K_Boxer_Tasha_1.0.bqsr.db.bed"
#
cd $FILE_DIR
#
find $PWD -maxdepth 2 -name "*L001*.bam" -printf '%f\n' &> "$tmpdir"/merge_files.tmp
find $PWD -maxdepth 2 -name "*L001*.bam" -printf '%h\n' &> "$tmpdir"/merge_dir.tmp
#
IFS=,$'\n' read -d '' -r -a lsdir < "$tmpdir"/merge_dir.tmp
#
for ((i = 0; i < ${#lsdir[@]}; i++))
do
	cd "${lsdir[$i]}";find . -name "*.bam" -printf '%f\n' | sort -V | paste -sd " " - >> "$tmpdir"/merge_lsfiles.tmp
done
#
cd $tmpdir
awk -F '_' '{print $1"_"$2"_"$3}' merge_files.tmp >> merge_IDs.tmp
awk -F '_' '{print $2"_"$3}' merge_files.tmp >> merge_samplename.tmp
#
IFS=,$'\n' read -d '' -r -a lsfiles < "$tmpdir"/merge_lsfiles.tmp
IFS=,$'\n' read -d '' -r -a ID < "$tmpdir"/merge_IDs.tmp 
IFS=,$'\n' read -d '' -r -a sample < "$tmpdir"/merge_samplename.tmp
#echo ${lsdir[@]}
#echo ${lsfiles[@]}
#echo ${ID[@]}
#echo ${sample[@]}
#sleep 30
#
for ((i = 0; i < ${#lsdir[@]}; i++))
do
	echo "cd ${lsdir[$i]}; samtools merge /lscratch/\$SLURM_JOB_ID/"${ID[$i]}".bam ${lsfiles[$i]} && gatk MarkDuplicates I=/lscratch/\$SLURM_JOB_ID/"${ID[$i]}".bam O="${lsdir[$i]}"/dedup_"${sample[$i]}".bam M="${lsdir[$i]}"/"${sample[$i]}".metrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true TMP_DIR=/lscratch/\$SLURM_JOB_ID && samtools index dedup_"${sample[$i]}".bam" >> "$homedir"/CF4_mergeDedup.swarm
done
#
cd $homedir
#
head CF4_mergeDedup.swarm
#
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
#jobid1=$(swarm -f CF4_mergeDedup.swarm -g 36 -t 20 --gres=lscratch:350 --time 2-0 --module samtools,GATK/4.1.9.0 --logdir ~/job_outputs/samtools/"$SWARM_NAME"_merge --sbatch "--mail-type=ALL,TIME_LIMIT_90 --job-name "$SWARM_NAME"_Merge")
#echo $jobid1
#
for ((i = 0; i < ${#lsdir[@]}; i++))
do
        echo "cd "${lsdir[$i]}"; gatk --java-options \"-Xmx16g -XX:ParallelGCThreads=4\" BaseRecalibrator -bqsr-baq-gap-open-penalty 40.0 -R "$CanFam4" --tmp-dir /lscratch/\$SLURM_JOB_ID -I dedup_"${sample[$i]}".bam --known-sites "$knownsite" -L "$INT" -O "${sample[$i]}"_firstpass_recal.table; gatk --java-options \"-Xmx16g -XX:ParallelGCThreads=4\" ApplyBQSR -R "$CanFam4" --tmp-dir /lscratch/\$SLURM_JOB_ID -I dedup_"${sample[$i]}".bam -bqsr "${sample[$i]}"_firstpass_recal.table -O "${sample[$i]}"_BQSR.bam" >> CF4_BQSR.swarm
done
#
head CF4_BQSR.swarm
read -p "Does the formatting of the swarmfile appear correct? (yes or no) " promptA
while true; do
        case "$promptA" in
                [YyEeSs]* ) break ;;
                [NnOo]* ) echo "Verify that the inputs are correct and try again"; exit;;
                *) echo "Enter yes or no" ;;
        esac
done
#
echo "BQSR Swarm Job ID: "
#
#jobid2=$(swarm -f bqsr_pipelinetest.swarm -g 36 -t 10 --time 4-0 --gres=lscratch:320 --module GATK/4.1.9.0 --logdir ~/job_outputs/gatk/BaseRecalibrator/"$SWARM_NAME"_BQSR --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME"_BQSR" --dependency:afterok:"$jobid1"")
#echo $jobid2
#
for ((i = 0; i < ${#lsdir[@]}; i++))
do
	echo "cd "${lsdir[$i]}"; gatk --java-options \"-Xmx16g -XX:ParallelGCThreads=4\" BaseRecalibrator -bqsr-baq-gap-open-penalty 40.0 -R "$CanFam4" --tmp-dir /lscratch/\$SLURM_JOB_ID -I "${sample[$i]}"_BQSR.bam --known-sites "$knownsite" -L "$INT2" -O "${sample[$i]}"_secondpass_recal.table; gatk AnalyzeCovariates --tmp-dir /lscratch/\$SLURM_JOB_ID -before "${sample[$i]}"_firstpass_recal.table -after "${sample[$i]}"_secondpass_recal.table -plots "${sample[$i]}"_Covariates.pdf" >> CF4_BQSR2pass_Covariates.swarm
done
#
#jobid3=$(swarm -f CF4_BQSR2pass.swarm -g 32 -t 8 --time 3-0 --gres=lscratch:250 --module GATK/4.1.9.0 --logdir ~/job_outputs/gatk/BaseRecalibrator/"$SWARM_NAME"_AC --sbatch "--mail-type=ALL,TIME_LIMIT_80 --job-name "$SWARM_NAME_AC" --dependency:afterok:"$jobid2"")
#
