# Preparing the coverage matrix in the JIC-HPC cluster

Author: Ricardo H. Ramirez Gonzalez

## Requirements

### Software
1. ```Bash``` and ```slurm```. Not required for the analysis, but the scripts are designed for an HPC environment. 
2. ```R```, with the [Uauy-Lab/bio.tilling](https://github.com/Uauy-Lab/bio.tilling) package
3. ```bedtools-2.28.0``` to extract the coverage of each cram file
4. ```samtools-1.9``` 


### Inputs

1. The ```fasta``` file, with the corresponding ```fai``` index. Used ```161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta```
2. A ```BED``` file with each of the regions sorrounding the genes. This is the file ```IWGSCv1.0_UTR_ALL_2000bp.merged.bed```
3. A set of ```BED``` files with the windows to use with the pattern: ```./bed_files/IWGSCv1.0_UTR_ALL_2000bp.merged.${bp}bp.bed ```. Where ```bp``` corresponds to the 100, 150 and 200 bp window sizes. 


## Description
The raw data is stored in a *cold* storage at the HPC. 

### 1. Extract the coverage from the reference

The script  ```slurm_downliad_cram_and_get_coverage.sh``` submits a job array to download the CRAM files to a local storage run ```bedCov```. The results are stored on a folder for each window (100, 150 and 200bp). 

Then, the script ```mergeBeds.sh``` joins the the coverage of each window in a single table

### 2. Normalise the table. 

The script ```detect_deletions_simple.R``` Calculates the normalisation of the matrix. 
The normalization goes both ways, by sample and by window. This is based on the ```detact_deletions.R``` script in ```Uauy-Lab/bio.tilling```, but it saves results as they are calculated.
 This is because the memory usage of the full matrix on this dataset his in the border of what we can use in the cluster. 

As we have 3 window sizes, the script ```submit_deletions.sh``` is used to submit each window as an array. 

### Preapring another filter 

Create a ```.genome``` file:

```
gunzip -c Triticum_aestivum.IWGSC.41.parts.tidy.gff3.gz | grep "##sequence-region" | awk '{print $2"\t"$4}' > Triticum_aestivum.IWGSC.41.parts.tidy.genome
```


Extract BED file only with the genes: 

```
gunzip -c Triticum_aestivum.IWGSC.41.parts.tidy.gff3.gz | grep IWGSC | grep ID=gene | awk '{print $1"\t"$4"\t"$5}' > Triticum_aestivum.IWGSC.41.parts.tidy.genes.bed
```


Get the expanded bed file using ```bedtools slop``` 


```
bedtools slop -g Triticum_aestivum.IWGSC.41.parts.tidy.genome -b 2000 -i Triticum_aestivum.IWGSC.41.parts.tidy.genes.bed | bedtools merge -i - > Triticum_aestivum.IWGSC.41.parts.tidy.genes.2k.merged.bed
```





