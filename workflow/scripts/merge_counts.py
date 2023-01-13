import pandas as pd
import glob

# Get a list of all files that end in "counts.csv"
countfiles = glob.glob("*counts.csv")
reffile = 'hg19_protein-coding_genes.bed'

# Create an empty list to hold the DataFrames
dfs = []

# Loop through the files and read them into a DataFrame
for file in countfiles:
    df = pd.read_csv(file)
    df.columns = ['Hugo_Symbol', str('NSM_Count_'+file).replace('.InGenes.counts.csv', '')]
    dfs.append(df)


# Merge the DataFrames on the "common_column"
merged_df = pd.read_csv(reffile, sep = '\t')
merged_df.columns = ['Chrom', 'Start', 'End', 'Hugo_Symbol', 'Length', 'Strand', '# Exons', 'Locus_Type', 'GC', 'ATGC']

for df in dfs:
        merged_df = merged_df.merge(df, on='Hugo_Symbol')

merged_df.to_csv('Merged.NSM_counts.csv')
