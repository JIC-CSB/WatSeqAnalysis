import logging

import dtoolcore

# Set up some logging.
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Data structure to hold unique accessions.
set_of_accessions = set()

# Create a storage broker to be able to list all the datasets in the base uri.
base_uri = "ecs://pr-raw-watseq"
storage_broker = dtoolcore._get_storage_broker(base_uri, None)

# Loop over all the datasets in the base uri.
for uri in storage_broker.list_dataset_uris(base_uri, None):

    # Load the dataset from the uri.
    dataset = dtoolcore.DataSet.from_uri(uri)

    # Load the overlay with the information about accessions.
    try:
        accession = dataset.get_overlay("accession")
    except dtoolcore.DtoolCoreKeyError:
        logger.info("{} has no accession overlay".format(dataset.name))
        continue

    # Loop over all the items in the base uri.
    for item_id in dataset.identifiers:

        # Only fastq items have an accession associated with them.
        if accession[item_id] is not None:
            set_of_accessions.add(accession[item_id])

for accession in set_of_accessions:
    print(accession)
