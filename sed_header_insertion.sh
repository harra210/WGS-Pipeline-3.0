#!/bin/bash
homedir=$(pwd)
cp final.txt header_test.txt
sed -i '1i IDfl IDln SM R1 R2 D1 D2 REF IDX FQDIR MAPDIR LOGDIR MEDIR QUALDIR' header_test.txt
head header_test.txt
