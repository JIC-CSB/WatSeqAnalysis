# Preparing the coverage matrix in the JIC-HPC cluster

Author: Ricardo H. Ramirez Gonzalez

## Requirements

### Software
1. ```Bash``` and ```slurm```. Not required for the analysis, but the scripts are designed for an HPC environment. 
2. ```R```, with the [Uauy/labbio.tilling](https://github.com/Uauy-Lab/bio.tilling) package
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





