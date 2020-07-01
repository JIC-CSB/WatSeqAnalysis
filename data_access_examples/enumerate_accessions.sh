# Location where all the watseq datasets are stored.
BASE_URI=ecs://pr-raw-watseq

# Loop over all the watseq datasets.
for URI in $(dtool ls -q $BASE_URI); do
    echo $URI;
    echo $(dtool name $URI)

    # Loop over all the items in a particular dataset.
    for ITEM_ID in $(dtool identifiers $URI); do

        # Only print out information about fastq files.
        IS_FASTQ=$(dtool item overlay is_fastq $URI $ITEM_ID)
        if [ ${IS_FASTQ,,} = "true" ]; then
            ACC=$(dtool item overlay accession $URI $ITEM_ID)
            RELPATH=$(dtool item relpath $URI $ITEM_ID)
            echo $ACC $ITEM_ID $RELPATH
        fi

    done
done
