#!/bin/bash
homedir=$(pwd)
cd tmp/
tmpdir=$(pwd)
files=( $(while IFS= read -r line; do head -n 1 "$line"; done < filldown.txt))
declare -a files
#
#echo ${files[@]}
#
cd $homedir
> files_test.txt
for f in ${files[@]}
do
	printf "%s " "$f" >> files_test.txt
done
cpRow=$(cat paste_test_add.txt | wc -l)
echo $cpRow
#more files_test.txt
#printf "%s\n" ${files[@]} >> paste_test_add.txt
> final.txt
for x in $cpRow
do
	awk '{for (i=1; i<=n; i++) print}' n="$x" files_test.txt &> files_test_awk.txt
done
more files_test_awk.txt
#sleep 10
#
paste paste_test_add.txt files_test_awk.txt > final.txt
sed -i $'s/\t/ /g' final.txt
#sed -i $'s/ /,/g' final.txt
more final.txt

