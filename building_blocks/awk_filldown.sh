#!/bin/bash
#using final.txt
> Final_Final.txt
cp final.txt final_awk_test.txt
#awk -F ',' -v OFS=',' '{ for (i = 1; i <= NF; i++) { if ($i == "") $i = save[i]; else save[i] = $i } for (i = NF+1; i <= 14; i++) save[i] = ""; print }' final_awk_test.txt &> Final_Final.txt
#awk -F',' -v OFS=',' \
#   '{
#        for (i = 1; i <= NF; i++) {
#                if ($i == "") $i = save[i]
#                else          save[i] = $i
#        }
#       for (i = NF+1; i <= 14; i++) save[i] = ""
#        print
#    }' final_awk_test.txt &> Final_Final.txt
#$@
#
#awk -F ',' -v OFS=',' 'FNR>1{if ($8!="") last=$8; else $8=last}1' final_awk_test.txt &> Final_Final.txt
#awk 'BEGIN{i=0;}{if ($2 == "" ){print i,$1;}else{print $0;}i=$1;}' final_awk_test.txt &> Final_Final.txt
#
#
#
#

more Final_Final.txt
