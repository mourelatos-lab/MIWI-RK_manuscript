#!/usr/bin/env python3
#
# Deduplicate tab separated file (two columns - first column is id, second column is vector)
# Keeps names of all duplicated rows (separated by |)
#

import sys
from collections import defaultdict

dedup_records = defaultdict(list)

if sys.argv[1] == "stdin":    
    ifile = sys.stdin
else:
    ifile = open(sys.argv[1], "r")

# separator for merging duplicated rows
if len(sys.argv) == 3:
    sep = sys.argv[2]
else:
    sep = "|"

for record in ifile:    
    record = record.rstrip().split('\t')
    # Use the sequence as the key and then have a list of id's as the value
    dedup_records[str(record[1])].append(record[0])
for seq, ids in dedup_records.items():
    sys.stdout.write(format(sep.join(ids)) + "\t" + seq + "\n")
