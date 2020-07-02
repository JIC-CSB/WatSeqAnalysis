#!/bin/bash
#SBATCH -p RG-Cristobal-Uauy,jic-medium,nbi-medium
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem 16000
#SBATCH -o log/extract_cov_\%A_\%a.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ricardo.ramirez-gonzalez@jic.ac.uk
#SBATCH --array=130
#SBATCH --time=2-00:00:00 
#SBATCH -J extract_cov
#SBATCH --localscratch=ssd:200

##   SBATCH --array=1- 8 24


pwd
hostname
date

source bedtools-2.28.0
source samtools-1.9

function log_line() {
	echo $(date) "$1" >&2
}

#this is the SSD folder. We requested 100 GB with --localscratch==ssd:100
SSD=$SLURM_LOCAL_SCRATCH

#use the ssd as dtool cache. This step is important because it avoids having the data on isilon
export DTOOL_CACHE_DIRECTORY=$SSD/dtool

mkdir -p $DTOOL_CACHE_DIRECTORY
echo $SSD
i=$SLURM_ARRAY_TASK_ID

let current_line=0

ref=../mappings/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta

tmp_d="$SSD/cram"
mkdir -p $tmp_d

chr3AFolder="../mappings/chr3A"
gene_bams="../mappings/around_genes"

mkdir -p "$chr3AFolder"
mkdir -p $gene_bams
cat metadata/crams_table.tsv | while read -r ecs prefix	cram_name cram_hash cram_index cram_index_hash; do

  if [[ $i -eq $current_line ]]; then

  	log_line "Loading $cram_name"
  	FILE=`dtool item fetch $ecs $cram_hash` 
  	mv $FILE $tmp_d/$cram_name
  	log_line "Finished $cram_name"


  	log_line "Loading $cram_index"
  	FILE=`dtool item fetch $ecs $cram_index_hash` 
  	mv $FILE $tmp_d/$cram_index
  	log_line "Finished $cram_index"

  	##For amber: 
  	log_line "Extracting: samtools view -b -T $ref -o $chr3AFolder/$prefix.bam $tmp_d/$cram_name chr3A_part1:26035170-26238727"
  	samtools view -b -T $ref -o $chr3AFolder/$prefix.bam $tmp_d/$cram_name chr3A_part1:26035170-26238727

  	##For the lab
  	log_line "Extracting: samtools view -C -T $ref -L ./bed_files/IWGSCv1.0_UTR_ALL_2000bp.merged.bed -o $gene_bams/$prefix.cram $tmp_d/$cram_name"
  	new_cram=$gene_bams/$prefix.cram
  	samtools view -C -T $ref -L ./bed_files/IWGSCv1.0_UTR_ALL_2000bp.merged.bed -o $new_cram $tmp_d/$cram_name 
  	

  	ls -lah $tmp_d
  	

  	for bp in 100 150 200; do 
  		./scripts/extractCoverageCRAM.sh $new_cram ./bed_files/IWGSCv1.0_UTR_ALL_2000bp.merged.${bp}bp.bed ./bedCov/window_${bp}bp $ref
  	done



  fi
  let current_line+=1
done 
du -hs $SSD/*
log_line "DONE"

#FILE=`dtool item fetch ecs://pr-raw-watseq/252a31bb-0cee-4359-afec-e8c87ed5810d 95fa22e416959144ea4d8b4835de9ebb6b4eaf6e` 
#mv $FILE ./merge.rmdup.B1190023-1.cram.crai

#FILE=`dtool item fetch ecs://pr-raw-watseq/252a31bb-0cee-4359-afec-e8c87ed5810d 38b705b4cd2b566c3f68410bd4617285f6be1790`
#mv $FILE ./merge.rmdup.B1190023-1.cram


#test of cram
#./scripts/extractCoverageCRAM.sh ../mappings/merge.rmdup.B1190023-1.cram ./bed_files/Triticum_aestivum.IWGSC.41_0bp_.bed ./bedCov ../mappings/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta

date
