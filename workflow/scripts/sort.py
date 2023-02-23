import pandas as pd
import sys

# read in the master table
MasterTable = pd.read_csv(sys.argv[1])

# sort the same way as previous table
MasterTable = MasterTable.sort_values(by=['Dist'])

# output the final table
MasterTable.to_csv(sys.argv[2], index=False, sep = '\t')