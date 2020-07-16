#!/bin/bash -e
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo $DIR
sqlite3 covs_200bp.db < $DIR/create_database.sql
gunzip -c lines.tsv.gz  | sqlite3   -csv -separator "	" covs_200bp.db '.import /dev/stdin lines'
gunzip -c windows.bed.gz  | sqlite3 -csv -separator "	" covs_200bp.db '.import /dev/stdin regions'
gunzip -c norm_covs.tsv.gz  | sqlite3 -csv -separator "	" covs_200bp.db '.import /dev/stdin norm_covs'
gunzip -c covs.tsv.gz  | sqlite3 -csv -separator "	" covs_200bp.db '.import /dev/stdin covs'