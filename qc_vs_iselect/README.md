# qc vs iselect

The purpose of this program is to check whole genome shotgun rawdata against genotyping data of a SNP chip. 

The basic idea is to convert a SNP into a set of kmers (in the manifest of the SNP chip you find flanking sequence for that) and then look for presence/absence of those kmers in the sequence data.

## input data

### Axiom® 35k Breeders' array probe set.

This is downloaded from [cerealsdb](https://www.cerealsdb.uk.net/cerealgenomics/CerealsDB/Excel/35K_array/35k_probe_set_IWGSCv1.xlsx) 
and converted to a tab separated text file.

I put a copy of it `35k_probe_set_IWGSCv1-1.txt` in this repository. 

### Genotypes of accessions

I got this file from Luzie Wingen. Not sure about copy restrictions.
The file `Watkins_829acc_axiom35k_gts_hmp.txt.gz` is in this repository. Unzip before use!

### A kmer table

A text file containing kmers from wgs sequencing. If a line in this file contains more than the kmer, it is considered to be tab separated and the kmer is the first column.

Generate such a file, for example, with jellyfish

```
source jellyfish-2.1.4
zcat accession*.fastq.gz | jellyfish count -t 20 -C -m 31 -s 4G -o accession.jf /dev/fd/0
jellyfish dump -L 3 -ct -o accession.dump.txt accession.jf
``` 


## The program

This is programmed with Java. Tested with jre-7.21. The code is under `src/`.
A runnable version is also available: `WatSeq_AxiomCheck.jar`.

Run with

```
mkdir output_dir
source jre-7.21
java -XX:ParallelGCThreads=2 -Xmx150000M -jar WatSeq_AxiomCheck.jar -m 35k_probe_set_IWGSCv1-1.txt -a Watkins_829acc_axiom35k_gts_hmp.txt -i accession.dump.txt -o output_dir
```

