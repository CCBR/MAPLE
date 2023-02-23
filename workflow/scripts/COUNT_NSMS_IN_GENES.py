# COUNT_NSMS_IN_GENES.py

import sys
import pandas as pd

def Count_NSMs_in_Genes(infile,sample,outfile):

    Hugo_Symbols = {}

    for line in open(infile, 'rt'):
        array = line.strip().split('\t')

        Hugo_Symbol = array[9]
        #Hugo_Symbol = array[7]

        if Hugo_Symbol not in Hugo_Symbols:
            Hugo_Symbols[Hugo_Symbol] = 0

        Hugo_Symbols[Hugo_Symbol] += 1


    Hugo_Symbols_df = pd.DataFrame.from_dict(Hugo_Symbols,orient='index')
    Hugo_Symbols_df.columns = [sample]


    Hugo_Symbols_df.to_csv(outfile)



Count_NSMs_in_Genes(sys.argv[1], sys.argv[2], sys.argv[3])
