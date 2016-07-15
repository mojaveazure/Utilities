#!/bin/bash

set -e
set -u
set -o pipefail

#   Define the usage message
function usage() {
    echo -e "\
Usage:  ./fakeHeader.sh <option> <file> \n\
\n\
where:  <option> is 'solo', 'serial', or 'parallel'  (without quotes) \n\
and:    <file> is a SAM or BAM file or text file listing SAM or BAM files \n\
\n\
This script adds the '@HD' line for a SAM file \n\
in case it was not added at the time of creation \n\
\n\
solo        Use 'solo' for adding the header line to one SAM file \n\
        When using 'solo', have <file> be a single SAM or BAM file \n\
\n\
serial      Use 'serial' to add the header line to multiple SAM files in serial \n\
        When using 'serial', have <file> be a text (.txt) file listing the full path to SAM files \n\
\n\
parallel    Use 'parallel' to add the header line to multiple SAM files in parallel \n\
        When using 'parallel', have <file> be a text (.txt) file listing the full path to SAM files \n\
        NOTE: This requires GNU Parallel to run \
" >&2
exit 1
}

#   Define a function to add the header line to a SAM file
function append_header() {
    sample_name=`basename "$1" .sam`
    sample_dir=`dirname "$1"`
    echo -e @HD"\t"VN:1.5"\t"SO:coordinate | cat - "$1" > "$sample_dir"/"$sample_name"_WithHeader.sam
}

#   Export the function for adding a header line
export -f append_header

#   Define a function to add the header line to a BAM file
function head_the_bam() {
    if `command -v samtools > /dev/null 2> /dev/null`
    then
        sample="$1"
        sample_name=`basename "$1" .bam`
        sample_dir=`dirname "$1"`
        samtools view -h "$sample" > "$sample_dir"/"$sample_name".sam
        echo -e @HD"\t"VN:1.5"\t"SO:coordinate | cat - "$sample_dir"/"$sample_name".sam > "$sample_dir"/"$sample_name"_WithHeader.sam
        samtools view -bhS "$sample_dir"/"$sample_name"_WithHeader.sam > "$sample_dir"/"$sample_name"_WithHeader.bam
        rm -f "$sample_dir"/"$sample_name".sam "$sample_dir"/"$sample_name"_WithHeader.sam
    else
        echo "SAMTools is required when using a BAM file"
        echo "Please install SAMTools and add to your PATH"
        exit 1
    fi
}

#   Export the function for adding a header line
export -f head_the_bam

#   Check to see if there are any arguments
#   If not, or less than the required 2, display the usage message
if [ "$#" -lt 2 ]
then
    usage
fi

#   Define arguments
METHOD="$1"
FILE="$2"

#   Check extension of input file, or files described in list, to determine function used
if [ "${FILE: -4}" == ".sam" ]
then
    FUNC=append_header
elif [ "${FILE: -4}" == ".bam" ]
then
    FUNC=head_the_bam
elif [ "${FILE: -4}" == ".txt" ]
then
    if grep -q "\.sam" "$FILE"
    then
        FUNC=append_header
    elif grep -q "\.bam" "$FILE"
    then
        FUNC=head_the_bam
    else
        echo "Incorrect file format in list"
        exit 1
    fi
else
    echo "Incorrect file format"
    exit 1
fi

#   Run the function for adding a header line according to the method
case "$METHOD" in
    "solo" )
        #   Only one file, run directly
        "$FUNC" "$FILE"
        ;;
    "parallel" )
        #   Run on multiple files in parallel
        #   Check to see if parallel is installed
        if `command -v parallel > /dev/null 2> /dev/null`
        then
            cat "$FILE" | parallel -v "$FUNC" {}
        else
            echo "GNU Parallel is not installed"
            echo "Please install or use 'serial' option'"
            exit 1
        fi
        ;;
    "serial" )
        #   Run on multiple files in serial (one after another)
        for i in `seq $(wc -l < "$FILE")`
        do
            sample=`head -"$i" "$FILE" | tail -1`
            "$FUNC" "$sample"
        done
        ;;
    * )
        #   For invalid methods, display usage information
        usage
        ;;
esac
