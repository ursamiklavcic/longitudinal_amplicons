#conda activate learn-python

import pandas as pd
impoert numpy as np

# Load CSV OTUTAB 
df = pd.read_csv(' G:\projekti\longitudinal_amplicons\dataqcsv_files\otutab_filt.csv')
otutab = df.to_numpy(df)

print(otutab)
