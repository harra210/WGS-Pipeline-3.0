#!/bin/bash
cd tmp
tmpdir=$(pwd)
cd $tmpdir
cd filldown
fill=$(pwd)
#
#files=( $(find $PWD -name "*.var" | sort -V | while IFS= read -r -a '' line; do; echo $ ))
find $PWD -name "*.var" | sort -V &> "$tmpdir"/filldown.txt
cd $fill
#mapfile -t files < filldown.txt
find . -name "*.var" -printf '%f\n' | sort -V |
	while IFS= read -r -d '' line
	do
		echo "$line"
	done
#IFS=,$'\n' read -d '' -r -a files < filldown.txt
#declare -a files
#echo ${files[@]}
#
#while IFS= read -r line;
#do
#	process "$line"
#done < filldown.txt
#cd $fill
#while read ${files[@]}
#do
#       sed -i "s/$/\t$f/" "$homedir"/paste_test_add.txt
#done < ${files[@]}
