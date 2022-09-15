# WeigthedDYADposition.py
# Find center(s) of fragments (DYADs) and give them weight (0.5 or 1)

import sys

infile = open(sys.argv[1], 'rt')
outfile = open(sys.argv[2], 'w+')

DYADs = {}

for line in infile:
    array = line.strip().split('\t')
    chrom, start, end = array[0:3]

    length = int(end) - int(start)

    # Length is a pair number. Center split in two positions

    if (length % 2) == 0:

        DYAD1 = int(int(start) + float(length)/2)
        DYAD2 = int(int(start) + float(length)/2 + 1)

        key1 = chrom+'|'+str(DYAD1)
        key2 = chrom+'|'+str(DYAD2)

        if key1 not in DYADs:
            DYADs[key1] = 0

        if key2 not in DYADs:
            DYADs[key2] = 0

        DYADs[key1] += 0.5
        DYADs[key2] += 0.5

    # Lenght is an odd number. Center is one position

    else:
        DYAD = int(int(start) + float(length)/2 + 0.5)

        key = chrom+'|'+str(DYAD)

        if key not in DYADs:
            DYADs[key] = 0

        DYADs[key] += 1

# Create output

for DYAD in DYADs:
    chrom, pos = DYAD.split('|')
    outfile.write(chrom+'\t'+pos+'\t'+str(DYADs[DYAD])+'\n')