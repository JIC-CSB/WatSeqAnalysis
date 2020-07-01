import logging

import dtoolcore
from dtool_info.utils import sizeof_fmt

# Set up some logging.
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Data for an accession can be spread over multiple dataset.
# Data structure to hold dataset URIs that have data from an accessions.
accession_dataset_uri_lookup = dict()

# Create a storage broker to be able to list all the datasets in the base uri.
base_uri = "ecs://pr-raw-watseq"
storage_broker = dtoolcore._get_storage_broker(base_uri, None)

# Loop over all the datasets in the base uri.
total_items = 0
total_bytes = 0
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
    num_items = 0
    num_bytes = 0
    for item_id in dataset.identifiers:

        acc = accession[item_id]
        # Only fastq items have an accession associated with them.
        if str(acc) != "None":
            num_items += 1
            num_bytes += dataset.item_properties(item_id)["size_in_bytes"]

    print("{}\t{}\t{}".format(dataset.name, num_items, sizeof_fmt(num_bytes)))
    total_items += num_items
    total_bytes += num_bytes

# Print out the summary statistics.
print("Total\t{}\t{}".format(total_items, sizeof_fmt(total_bytes)))
