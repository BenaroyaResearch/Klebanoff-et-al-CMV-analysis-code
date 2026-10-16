

### contig processing for sam/ campbell lab samples 

## Goal - create a compiled product with barcode-tagged TRA/TRB pairs and queue to timepoint 


## to run on the terminal - ### not running this on gizmo was a bad idea. 
## ssh rhino02
## source ~/.zshrc
## conda activate file_mover
# cd /fh/fast/corey_l/esford3_kmayerbl_collab/cmv_tcrs/klebanoff_2025/

ipython
import os
import numpy as np
import pandas as pd
from tenextra.parse import select_likely_receptors


## handle contigs
store_ac = []
store_cc = []

files = os.listdir()

for file in files: 
  if file.endswith('filtered_contig_annotations.csv'): 
    print(file)
    clean_clones, all_clones, ct_chains, ct_pairs = \
      select_likely_receptors(
          f = file,
          threshold_chains = 10)
    clean_clones['filename'] = file
    all_clones['filename'] = file
    store_cc.append(clean_clones)
    store_ac.append(all_clones)

# combine into a single dataframe
all_clean_clones = pd.concat(store_cc, ignore_index = True)
all_all_clones = pd.concat(store_ac, ignore_index = True)

all_clean_clones['v_a_gene']=all_clean_clones['v_a_gene'].apply(lambda x: f"{x}*01")
all_clean_clones['j_a_gene']=all_clean_clones['j_a_gene'].apply(lambda x: f"{x}*01")
all_clean_clones['v_b_gene']=all_clean_clones['v_b_gene'].apply(lambda x: f"{x}*01")
all_clean_clones['j_b_gene']=all_clean_clones['j_b_gene'].apply(lambda x: f"{x}*01")

all_clean_clones.columns

#In [12]: all_clean_clones.shape
#Out[12]: (75241, 69)

df = all_clean_clones[['barcode', 'v_a_gene','v_b_gene','j_a_gene','j_b_gene','cdr3_a_aa','cdr3_b_aa',
                                                'cdr3_a_nt','cdr3_b_nt','filename']]

df.to_csv('all_donor_barcode_tcrs_121325.csv', index = False)


