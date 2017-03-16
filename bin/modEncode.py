#!/usr/bin/env python

from intermine.webservice import Service
service = Service("http://intermine.modencode.org/release-33/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Submission")

# The view specifies the output columns
query.add_view(
    "DCCid", "description", "experimentType", "qualityControl",
    "databaseRecords.accession"
)

# This query's custom sort order is specified below:
query.add_sort_order("Submission.experimentType", "ASC")

# You can edit the constraint values below
query.add_constraint("databaseRecords.database", "=", "SRA", code = "A")
query.add_constraint("databaseRecords.accession", "CONTAINS", "SRR", code = "B")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A and B")

for row in query.rows():
    print(row["DCCid"], row["description"], row["experimentType"], row["qualityControl"], row["databaseRecords.accession"])
