#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 13:31:10 2021

@author: user
"""


import os
import glob
import pandas as pd
import numpy as np

file_path='/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/HIOE_NEUROG3_induction'

#==============================================
# add annotated ChIP onto TF-target table on specific timepoint
#==============================================
tfl=['NEUROG3','NEUROD1','RFX6']
cond_time=24
weight=1

rna_de = pd.read_table(f'/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/XZ_07-21-28-08-05-2020_RNA-seq_HIOE/Counts_TPM_mat/analysis/results/WT_HIOE_{cond_time}hpi_filtered.txt',sep='\t', index_col=0).index.tolist()

all_path=f'{file_path}/outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_x10_bias10_maxComb/{cond_time}hpi_Cores/Core_prior_atac_Miraldi_q_ChIP_x10_bias10_maxComb_fdr5_HIOE_NEUROG3_{cond_time}hpiSet_All.tsv'
de_path=f'{file_path}/outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_x10_bias10_maxComb/{cond_time}hpi_Cores/Core_prior_atac_Miraldi_q_ChIP_x10_bias10_maxComb_fdr5_HIOE_NEUROG3_{cond_time}hpiSet_DE.tsv'

all_df=pd.read_csv(all_path,sep='\t',index_col=0)
de_df=pd.read_csv(de_path,sep='\t',index_col=0)

for p, prior in enumerate([all_df, de_df]):
    prior_df=prior
    
    if p==0:
        gset='All_ChIP'
    else:
        gset='DE_ChIP'

    for tf in tfl:
        for f in glob.glob(f'{file_path}/inputs/priors/ChIP_annotated/{tf}*annota*.txt'):
            chip_anno=pd.read_table(f)
            gene_count={k:1 for k in set(chip_anno['Gene Name'].tolist())}

            if tf not in prior_df.columns:
                prior_df[f'{tf}']=[0]*prior_df.shape[0]
            
            for k in gene_count:
                if p==1 and not(k in rna_de):
                    continue
                elif k in prior_df.index:
                    prior_df.loc[k,f'{tf}']=(prior_df.loc[k,f'{tf}']+gene_count[k]*weight)
                    
    if weight!=1:
        prior_df.to_csv(f'{os.path.dirname(all_path)}/Core_prior_atac_Miraldi_q_ChIP_x10_bias10_maxComb_fdr5_HIOE_NEUROG3_{cond_time}hpiSet_{gset}_x{weight}.tsv',sep='\t',header=True)
    else:
        prior_df.to_csv(f'{os.path.dirname(all_path)}/Core_prior_atac_Miraldi_q_ChIP_x10_bias10_maxComb_fdr5_HIOE_NEUROG3_{cond_time}hpiSet_{gset}.tsv',sep='\t',header=True)


