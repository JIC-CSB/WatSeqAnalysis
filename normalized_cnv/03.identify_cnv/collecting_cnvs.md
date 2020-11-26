



# Reducing the number of columns to store teh database


```sh
gunzip -c all_cnv.csv.gz | awk -F"," '{print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8}' | gzip -c > all.cnv.for.db.tsv.gz 

```

The header for the table is: ```seqnames	start	end norm_cov	cnv_level	line```

seqnames,start,end,width,strand,norm_cov,cnv_level,line,max_gap


#To run the stiching program

```sh
 Rscript bin/stich_cnvs.r -i /Volumes/ExtremeSSD/watseq/window_coverage/bedCov/deletions_20200724/200bp/long_tables/cnv_out_M10/all_cnv.csv.gz -o /Volumes/ExtremeSSD/watseq/window_coverage/bedCov/deletions_20200724/200bp/long_tables/cnv_out_M10/M10_cnv  -g /Users/ramirezr/Dropbox/JIC/WatSeq/Annotation/IWGSCv1.0_UTR_ALL.parts.tidy.gff3.gz --cores 1 --chromosome chr6B_part2
```

