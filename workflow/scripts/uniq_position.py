# Uniq_Position.py
# Make histogram of DYADs

import sys

infile = open(sys.argv[1], 'rt')
outfile = open(sys.argv[2], 'w+')

last_key = 'chr1|1'
sumcount = 0

# Sum weights at each position

for line in infile:
    chrom, position, count = line.strip().split('\t')

    key = chrom + '|' + position

    if key == last_key:
        sumcount += float(count)

    else:
        if sumcount > 0:
            last_chrom, last_position = last_key.split('|')
            outfile.write(last_chrom+'\t'+last_position+'\t'+str(sumcount)+'\n')
        last_key = key
        sumcount = float(count)

last_chrom, last_position = last_key.split('|')

# Create output

outfile.write(last_chrom+'\t'+last_position+'\t'+str(sumcount)+'\n')