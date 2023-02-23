import pandas as pd
import glob
import sys
from os import listdir
from os.path import isfile, join

# Get a list of all files that end in "counts.csv"
countfiles = [f for f in listdir(sys.argv[1]) if isfile(join(sys.argv[1], f))]
reffile = sys.argv[2]
outputfile = sys.argv[3]

# Create an empty list to hold the DataFrames
dfs = []

# Loop through the files and read them into a DataFrame
for file in countfiles:
    df = pd.read_csv(sys.argv[1]+ "/" + file)
    df.columns = ['Hugo_Symbol', str('NSM_Count_'+file).replace('.InGenes.counts.csv', '')]
    dfs.append(df)

# Merge the DataFrames on the "common_column"
merged_df = pd.read_csv(reffile, sep = '\t')
merged_df.columns = ['Chrom', 'Start', 'End', 'Hugo_Symbol', 'Length', 'Strand', '# Exons', 'Locus_Type', 'GC', 'ATGC']

for df in dfs:
        merged_df = merged_df.merge(df, on='Hugo_Symbol', how="outer")

# sort df by SYMBOL
sorted_df=merged_df.sort_values(by=['Hugo_Symbol'], ascending=True)
sorted_df.reset_index(drop=True, inplace=True)

# output file
sorted_df.to_csv(outputfile)
