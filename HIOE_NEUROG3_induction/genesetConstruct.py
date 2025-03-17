#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 15:05:37 2020

@author: user
"""

import os
import pandas as pd

os.chdir('/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/HIOE_NEUROG3_induction')

# Construct geneset for each condition
for t in [24,48,72,96]:
    gene_set_file=f'/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/XZ_07-21-28-08-05-2020_RNA-seq_HIOE/Counts_TPM_mat/analysis/results/WT_HIOE_{t}hpi_filtered.txt'
    # geneset output FILE NAME cannot contain _
    geneset_condition=f'{t}hpiSet'
    geneset_up_id=f'{t}hpi_up'
    geneset_down_id=f'{t}hpi_down'
    
    gene_set=pd.read_table(gene_set_file,sep='\t')
    geneset_up=gene_set[gene_set.log2FoldChange>0].Gene.tolist()
    geneset_down=gene_set[gene_set.log2FoldChange<0].Gene.tolist()
    
    geneset_df=pd.DataFrame(data={'set_id':[geneset_up_id,geneset_down_id],
                                  'set_name':[geneset_up_id,geneset_down_id],
                                  'gene_set':['|'.join(geneset_up),'|'.join(geneset_down)]})
    
    geneset_df.to_csv(f'./inputs/geneSets/{geneset_condition}.txt',sep=' ', index=False, header=False)

# Construct leaveout list
gene_exp='inputs/geneExpression/RNAseq_24_DESeq2_VSDcounts.txt'
gene_exp_list=pd.read_table(gene_exp,sep='\t', index_col=0).columns.tolist()

for cond in [0,100]:
    for t in [24,48,72,96]:
        leveout=[x for x in gene_exp_list if ('_'+str(cond)+'_' in x) and (str(t)+'hpi' in x)]
        with open (f'inputs/leaveOutLists/{cond}-{t}hpi.txt','w') as f:
            f.write('\n'.join(leveout))