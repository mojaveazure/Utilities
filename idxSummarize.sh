#!/bin/bash

set -e
set -o pipefail

#   Check to make sure that datamash and parallel are installed
if ! $(command -v datamash > /dev/null 2> /dev/null); then echo "Please install GNU datamash" >&2; exit 1; fi
if ! $(command -v parallel > /dev/null 2> /dev/null); then echo "Please install GNU Parallel" >&2; exit 1; fi

#   Create a usage message
function Usage() {
    echo -e "\
Usage:  $(basename $0) <stats_list> [project] \n\
Where:  <sample_list>   is a list of SAMTools idxstats output files\n\
        [project]       is an optional name for the output file, defaults to 'STATS' \n\
" >&2
    exit 1
}

#   Export the function
export -f Usage

#   A function to sum the sequence length, number of mapped reads, and number of unmapped reads for a single sample
function getCounts() {
    local sample="$1" # What sample are we working on?
    local sampleExtension=".$(echo ${sample} | rev | cut -f 1 -d '.' | rev)" # Get the extension
    local sampleName="$(basename ${sample} ${sampleExtension})" # Remove the extension and file path
    local sequenceLength=$(cat ${sample} | head -$(($(wc -l < ${sample})-1)) | datamash sum 2) # Get the total sequence length
    local mappedReads=$(cat ${sample} | head -$(($(wc -l < ${sample})-1)) | datamash sum 3) # Get the total number of mapped reads
    local unmappedReads=$(cat ${sample} | head -$(($(wc -l < ${sample})-1)) | datamash sum 4) # Get the total number of unmapped reads
    echo -e "${sampleName}\t${sequenceLength}\t${mappedReads}\t${unmappedReads}" # Return the values, tab-delimeted
}

#   Export the function
export -f getCounts

#   Ensure we have our required one argument
if [[ "$#" -lt 1 ]]; then Usage; fi

#   Collect the arguments
STATS_LIST="$1"
PROJECT="${2:-'STATS'}"

#   Check to make sure our stats list exists
if ! [[ -f "${STATS_LIST}" ]]; then echo "Cannot find ${STATS_LIST}" >&2; exit 1; fi

#   Get the stats directory
STATS_DIR=$(dirname "${STATS_LIST}")

#   Create a name for our output
OUTPUT="${STATS_DIR}/${PROJECT}_allStats.txt"

#   Run getCounts in parallel and write to output
parallel getCounts {} > "${OUTPUT}" :::: "${STATS_LIST}"
