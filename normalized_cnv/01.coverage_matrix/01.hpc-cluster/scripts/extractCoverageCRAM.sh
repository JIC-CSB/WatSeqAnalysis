#!/bin/bash -e

inputCRAM=$1
exonsFile=$2
outputFolder=$3
fasta=$4



function log_line() {
	echo $(date) "$1" >&2
}

fai=$fasta.fai

source bedtools-2.28.0
source samtools-1.9
f=$inputCRAM
mkdir -p $outputFolder

filename=$(basename -- "$inputCRAM")
extension="${filename##*.}"
filename="${filename%.*}"

out_bedCov="$outputFolder/$filename.bedCov.gz"
if [ -f $out_bedCov ] ; then
    log_line "Already extracted ($out_bedCov)"
else
	log_line "running CovBed for $exonsFile $f and $fai to save in $out_bedCov"
    samtools view -b -T $fasta $inputCRAM | bedtools coverage -a $exonsFile -b stdin -sorted -g $fai | gzip -c  > $out_bedCov 
    log_line "Finished extracting $out_bedCov"
fi

