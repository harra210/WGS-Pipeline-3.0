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
# Colourise the output
RED='\033[0;31m'        # Red
GRE='\033[0;32m'        # Green
YEL='\033[1;33m'        # Yellow
NCL='\033[0m'           # No Color

file_specification() {
        FILE_NAME="$(basename "${entry}")"
        DIR="$(dirname "${entry}")"
        NAME="${FILE_NAME%.*}"
        EXT="${FILE_NAME##*.}"
        SIZE="$(du -sh "${entry}" | cut -f1)"

        printf "%*s${GRE}%s${NCL}\n"                    $((indent+4)) '' "${entry}"
        printf "%*s\tFile name:\t${YEL}%s${NCL}\n"      $((indent+4)) '' "$FILE_NAME"
        printf "%*s\tDirectory:\t${YEL}%s${NCL}\n"      $((indent+4)) '' "$DIR"
        printf "%*s\tName only:\t${YEL}%s${NCL}\n"      $((indent+4)) '' "$NAME"
        printf "%*s\tExtension:\t${YEL}%s${NCL}\n"      $((indent+4)) '' "$EXT"
        printf "%*s\tFile size:\t${YEL}%s${NCL}\n"      $((indent+4)) '' "$SIZE"
}

walk() {
        local indent="${2:-0}"
        printf "\n%*s${RED}%s${NCL}\n\n" "$indent" '' "$1"
        # If the entry is a file do some operations
        for entry in "$1"/*; do [[ -f "$entry" ]] && file_specification; done
        # If the entry is a directory call walk() == create recursion
        for entry in "$1"/*; do [[ -d "$entry" ]] && walk "$entry" $((indent+4)); done
}

# If the path is empty use the current, otherwise convert relative to absolute; Exec walk()
[[ -z "${1}" ]] && ABS_PATH="${PWD}" || cd "${1}" && ABS_PATH="${PWD}"
walk "${ABS_PATH}"      
echo     
