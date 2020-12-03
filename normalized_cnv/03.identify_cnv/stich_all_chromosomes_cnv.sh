#!/bin/bash

project_folder=/Volumes/ExtremeSSD/watseq/window_coverage/bedCov/deletions_20200724
m=M10
cnvs_file=$project_folder/200bp/long_tables/cnv_out_$m/all_cnv.csv.gz
gff=/Users/ramirezr/Dropbox/JIC/WatSeq/Annotation/IWGSCv1.0_UTR_ALL.parts.tidy.gff3.gz
out=$project_folder/200bp/long_tables/cnv_out_$m/${m}_cnv
cores=1
biocnv=/Users/ramirezr/Documents/public_code/WatSeqAnalysis/normalized_cnv/bio.cnv
cd $biocnv
for f in `gunzip -c $gff | grep sequence | head -42 | awk '{print $2}'`; do
#for f in `gunzip -c $gff | grep sequence | head -1 | awk '{print $2}'`; do 
    echo $f
     Rscript ./bin/stich_cnvs.r -i $cnvs_file -o $out  -g $gff --cores $cores --chromosome $f
done