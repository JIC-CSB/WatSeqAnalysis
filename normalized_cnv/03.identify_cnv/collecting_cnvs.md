



# Reducing the number of columns to store teh database


```sh
gunzip -c all_cnv.csv.gz | awk -F"," '{print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8}' | gzip -c > all.cnv.for.db.tsv.gz 

```

The header for the table is: ```seqnames	start	end norm_cov	cnv_level	line```

seqnames,start,end,width,strand,norm_cov,cnv_level,line,max_gap

