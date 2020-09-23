#!/bin/bash
> sample_swarmfile.txt
homedir=$(pwd)
while IFS=" " read -r FC TAG LB
do
	echo "cd "$homedir"; RGID=10701 RGPU="$FC"."$TAG" RGSM=10701 RGPL=ILLUMINA RGLB="$LB"" >> "$homedir"/sample_swarmfile.txt
done < pasted_header.txt
