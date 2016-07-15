#!/bin/bash

set -e
set -u
set -o pipefail

#   A script to subsample Fasta/FastQ files using Seqtk
if ! `command -v seqtk > /dev/null 2> /dev/null`
then
    echo "Failed to find seqtk!"
    exit 1
fi

#   Define a usage message
function Usage() {
    echo -e "\
Usage:  subsample.sh sample_info frac [outdir] [seed]\n\
\n\
where:  'sample_info' is a Fasta/FastQ file to be subsampled \n\
            or a list of Fasta/FastQ files to be subsampled \n\
        'frac' is either the fraction or number of reads to be kept \n\
\n\
Optional:   'outdir' is the directory to put the subsampled files \n\
                If not specified, we will use \n\
                `pwd`
\n\
            'seed' is the seed value to be used for all samples run during \n\
                this instance. If not specified, a random number will \n\
                be generated and recorded for posterity's sake \n\
" >&2
    exit 1
}

#   Do we have the two required arguments?
if [ "$#" -lt 2 ]
then
    Usage
fi

#   Give the arguments names and set default values
SAMPLE_INFO=$1
FRAC=$2
OUTDIR=${3:-`pwd`}
SEED=${4:-`echo "$RANDOM"`}

#   Make sure that the outdirectory exists
mkdir -p ${OUTDIR}

#   Define a function to carry out the subsampling
function subsample() {
    sample="$1"
    seed="$2"
    frac="$3"
    out="$4"
    echo "Subsampling ${sample}"
    name=`basename "${sample}" | cut -d '.' -f 1`
    ext=`basename "${sample}" | cut -d '.' -f 2`
    seqtk sample -s"${seed}" "${sample}" "${frac}" > "${out}"/"${name}"_sub_"${frac}"."${ext}"
}

#   Export this function
export -f subsample

#   Get the first line of ${SAMPLE_INFO}
LINE=`head -1 "${SAMPLE_INFO}"`

#   See if this first line is a file
if [[ -f "${LINE}" ]]
then #  If so
    if `command -v parallel > /dev/null 2> /dev/null` # See if we have GNU Parallel installed
    then #  If so
        cat "${SAMPLE_INFO}" | parallel "subsample {} ${SEED} ${FRAC} ${OUTDIR}" #  Run the subsample function on all files in parallel
    else #  If we don't have parallel
        for fastq in `cat "${SAMPLE_INFO}"` #   Run the subsample function in serial
        do
            subsample "$fastq" "${SEED}" "${FRAC}" "${OUTDIR}"
        done
    fi
else #  If the first line doesn't point to a file, assume we were given a Fasta/FastQ file
    subsample "${SAMPLE_INFO}" "${SEED}" "${FRAC}" "${OUTDIR}"
fi

#   Record what ${SEED} was, let us know where the final files can be found, and write all this to an out message file
echo -e "`date`:\nWe used a seed value of ${SEED}\nAll subsampled files can be found at ${OUTDIR}\nThis message can be read again at ${OUTDIR}/outMSG.txt" | tee -a ${OUTDIR}/outMSG.txt
