# Usage:  get_fastq_paired_reads.sh accession dataset_uri [dataset_uri]
# The first argument to the script is the accession.
# The remianing arguments to the script are dataset URIs where the accession is stored.

ACCESSION=$1
DATASET_URIS=${@:2}

# Make sure that the dtool cache is set to be in the /jic/scratch/projects/watseq directory.
export DTOOL_CACHE_DIRECTORY=/jic/scratch/projects/watseq/dtool_cache

for URI in $DATASET_URIS; do
    NAME=$(dtool name $URI)
    echo "Getting accessions from $NAME $URI";

    # Loop over all the items in a particular dataset.
    for ITEM_ID in $(dtool identifiers $URI); do

        # Only work on fastq files associated with the accession.
        ACC=$(dtool item overlay accession $URI $ITEM_ID)
        if [ $ACC = $ACCESSION ]; then

            # Only work on read1 fastq files.
            IS_READ1=$(dtool item overlay is_read1 $URI $ITEM_ID)
            if [ ${IS_READ1,,} = "true" ]; then

                # Get the item identifier of read2.
                PAIR_ID=$(dtool item overlay pair_id $URI $ITEM_ID)

                # Get the relative paths.
                READ1_RELPATH=$(dtool item relpath $URI $ITEM_ID)
                READ2_RELPATH=$(dtool item relpath $URI $PAIR_ID)

                
                # Do some logggin.
                echo $ACC $ITEM_ID $PAIR_ID $READ1_RELPATH $READ2_RELPATH

                # Get absolute paths to where the contents of the reads will be located.
                # Uncomment the lines below to get access to absolute paths
                # where the content of the read data will be available.
                #READ1_ABSPATH=$(dtool item fetch $URI $ITEM_ID)
                #READ2_ABSPATH=$(dtool item fetch $URI $ITEM_ID)
            fi

        fi

    done
done
