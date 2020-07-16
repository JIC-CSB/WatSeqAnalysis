# 02. Prepare the SQLite database

Author: Ricardo H. Ramirez Gonzalez

## Requirements

### Software
1. ```SQLite3``` To store the database
2. ```Ruby``` to run the
3. ``bash``

### Inputs

1. The files ```mat.csv.gz``` and ```df.csv.gz```, produced from the outputs in 01. This files where compressed with ```gzip```


## Description

To be able to do the analysis without keeping the full normalized coverage matrix in memory, I developed a simple ```sqlite``` database with three tables: 

1. ```lines``` The list of all the lines in the matrix
2. ```regions``` Each of the windows 
3.  ```norm_covs``` The melted values for each line and region. 


### 01. Prepare the tables to be imported in sqlite3

Because the matrix of coverage is large, the best approach is to avoid keeping everything in memory and to avoid lading the table on idividual statements. 
To take advantage of the fact that the ```mat.csv.gz``` and ```df.csv.gz``` come in the same order, the ruby script ```meltMat.rb``` takes a folder containing both of the scripts printes the ```norm_covs.tsv``` and ```windows.bed``` by melting each line in the matrix. Each line and window size get an ```INT``` id that is used to join the tables and to keep the files as small as possible.
The last file to be printed, is ```lines.tsv``` which prints the line names, based on the column name in the matrix. 
It is critical that this files are created with this script as it ensures the consistency of the primary keys. 

### 02. Import the database

```sh
sqlite3 covs_200bp.db < create_tables_sqlite.sql
cat lines.tsv  | sqlite3   -csv -separator "	" covs_200bp.db '.import /dev/stdin lines'
gunzip -c windows.bed.gz  | sqlite3 -csv -separator "	" covs_200bp.db '.import /dev/stdin regions'
gunzip -c norm_covs.tsv.gz  | sqlite3 -csv -separator "	" covs_200bp.db '.import /dev/stdin norm_covs'
gunzip -c covs.tsv.gz  | sqlite3 -csv -separator "	" covs_200bp.db '.import /dev/stdin covs'
```

