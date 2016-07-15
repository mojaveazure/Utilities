#!/bin/bash

set -e
set -u
set -o pipefail

#   A simple script to calculate coverage depth
#   from the text files from BEDTools coverage

function Usage() { #    Display a usage message and exit
    echo -e "\
Usage:  ./getCoverage.sh sample_info [outfile] \n\
where:      sample_info is a list of coverage histograms from BEDTools coverage \n\
\n\
Optionally: You can specify a name for the outfile \n\
    The default is 'allCoverage.txt' \n\
\n\
This script is dependent on GNU Parallel \n\
    Please install to use \n\
" >&2
    exit 1
}

#   Check to see if GNU Parallel is installed
if ! `command -v parallel > /dev/null 2> /dev/null`
then
    echo "GNU Parallel not found!"
    exit 1
fi

#   Make sure there's at least one argument passed to the script
#   otherwise, run the Usage function
if [ "$#" -lt 1 ]
then
    Usage
fi

#   Assign variable names to arguments
SAMPLE_INFO=$1
OUTFILE=${2:-allCoverage.txt}

#   Change to the directory where SAMPLE_INFO is
cd `dirname "${SAMPLE_INFO}"`

#   Write a function to calculate the coverage per sample
function coverageDepth() {
    SAMPLE="$1"
    NAME=`basename "${SAMPLE}" .coverage.hist.txt`
    echo "Searching for the depths"
    declare -ai depth=(`grep 'all' "${SAMPLE}" | cut -f 2`) #   Create an array containing the depths
    echo "Searching for the fraction of ${NAME} at each depth"
    declare -ai frac=(`grep 'all' "${SAMPLE}" | cut -f 5`) #    Create an array containing the fractions
    if ! [ ${#depth[@]} -eq ${#frac[@]} ] # Check to make sure the arrays are of equal length
    then
        echo "Failed to find equal amounts of coverage depth and fractions!"
        exit 1
    else
        echo "Calculating the coverage for ${NAME}"
    fi
    iterCov=0 # Start an object to calculate the entire coverage for the sample
    for i in `seq 0 $[ ${#depth[@]} - 1 ]`
    do
        #   Use 'bc' to allow for floating point calculations
        coverage=`echo "${depth[$i]} * ${frac[$i]}" | bc` # Calculate coverage per depth
        iterCov=`echo "${iterCov} + ${coverage}" | bc` #    Add to the previous cumalative coverage
    done
    echo -e "${NAME}:\t${iterCov}" > "${NAME}"_cov.txt #    Create a holding file for each sample
}

export -f coverageDepth #   Export the function to be used by GNU Parallel

cat ${SAMPLE_INFO} | parallel coverageDepth {} #    Run the coverageDepth function in parallel

find `pwd` -name "*_cov.txt" > allSamples.cov.txt # Find all the files made

cat allSamples.cov.txt | parallel cat {} > ${OUTFILE} #    Concatenate them into one file

rm -rf *cov.txt #Remove the intermediate files
