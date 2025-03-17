# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import re
import glob
import pandas as pd
import numpy as np

file_path='/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/HIOE_NEUROG3_induction/inputs'

#==============================================
# gtf to bed file as gene body
#==============================================
!gunzip -c gencode.v34lift37.annotation.gtf.gz |  awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$14,$7}}' | tr -d '";' | head
!gunzip -c gencode.v34lift37.annotation.gtf.gz |  awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$14,$7}}' | tr -d '";' > gene_body.bed

#==============================================
# merge ATAC peaks
#==============================================
#!cat *.bed | bedtools sort | bedtools merge > output.bed

#==============================================
# TFs to list, intersect with RNA-seq expression
# normGene: vsd normalized gene x sample matrix
# genesforTFA: ALL genes from RNA-seq across samples
# targetGene: Differentially expressed genes from RNA-seq across samples
# potRegs: Differentially expressed genes from RNA-seq across samples INTERSECT with motif2gene list
#==============================================
tfs_anno=pd.read_table(f'{file_path}/geneSets/human_TF_set.txt')
rna=pd.read_table(f'{file_path}/geneExpression/RNAseq_24_DESeq2_VSDcounts.txt')

targetGene=set()
for f in glob.glob('/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/XZ_07-21-28-08-05-2020_RNA-seq_HIOE/Counts_TPM_mat/analysis/results/*hpi_filtered.txt'):
    de_pd=pd.read_table(f)
    de_gene=de_pd.Gene.tolist()
    targetGene=targetGene | set(de_gene)
with open(f'{file_path}/targRegLists/targetGenes_names.txt','w') as f:
    f.writelines('\n'.join(targetGene))
    
potRegs=list(set(tfs_anno.Gene.tolist()) & targetGene)
tf_add='NEUROG3'
if tf_add not in potRegs:
    potRegs.append(tf_add)    
with open(f'{file_path}/targRegLists/potRegs_names.txt','w') as f:
    f.writelines('\n'.join(potRegs))

# background (DEFAULT all genes, optionally slice relieved padj as background) and genesForTFA
genesForTFA=rna
genesForTFA.to_csv(f'{file_path}/targRegLists/genesForTFA.txt',sep='\t',header=True,index=False)
genesForTFA.iloc[:,0].to_csv(f'{file_path}/geneSets/background_set.txt',sep='\t',header=False,index=False)

# background with DE
for p in [0.05, 0.15, 0.25]:
    bg_gene=set() 
    for f in glob.glob('/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/XZ_07-21-28-08-05-2020_RNA-seq_HIOE/Counts_TPM_mat/analysis/results/*hpi.txt'):
        de_pd=pd.read_table(f)
        de_gene=de_pd[(abs(de_pd.log2FoldChange)>1) & (de_pd.padj<p)].Gene.tolist()
        bg_gene=bg_gene | set(de_gene)
    with open(f'{file_path}/geneSets/background_set_de{p}.txt','w') as f:
        f.writelines('\n'.join(bg_gene))
#==============================================
# add annotated ChIP onto merged ATAC priors
#==============================================
tfl=['NEUROG3','NEUROD1','RFX6']
weight=1

prior_df_q=pd.read_csv(f'{file_path}/priors/prior_atac_Miraldi_q.tsv',sep='\t',index_col=0)
prior_df_b=pd.read_csv(f'{file_path}/priors/prior_atac_Miraldi_b.tsv',sep='\t',index_col=0)

for p, prior in enumerate([prior_df_q, prior_df_b]):
    prior_df=prior
    
    if p==0:
        for tf in tfl:
            for f in glob.glob(f'{file_path}/priors/ChIP_annotated/{tf}*annota*.txt'):
                chip_anno=pd.read_table(f)
                ps_mean=chip_anno.groupby('Gene Name')['Peak Score'].mean()
                if ps_mean.mean()!=0:
                    # weighted log(mean peak score) for each gene called by EXPERIMENTAL ChIP
                    genel=np.log2(ps_mean+1)
                else:
                    genel=chip_anno.groupby('Gene Name')['Peak Score'].count()
                gene_count=genel.to_dict()
    
                if tf not in prior_df.columns:
                    prior_df[f'{tf}']=[0]*prior_df.shape[0]
                
                for k in gene_count:
                    if k not in prior_df.index:
                        prior_df.append({f'{tf}':gene_count[k]*weight},ignore_index=True).fillna(0)
                    else:
                        prior_df.loc[k,f'{tf}']=prior_df.loc[k,f'{tf}']+gene_count[k]*weight
                        
        if weight!=1:
            prior_df.to_csv(f'{file_path}/priors/prior_atac_Miraldi_q_ChIP_x{weight}.tsv',sep='\t',header=True)
        else:
            prior_df.to_csv(f'{file_path}/priors/prior_atac_Miraldi_q_ChIP.tsv',sep='\t',header=True)
                
                                                     
    if p==1:
        for tf in tfl:
            for f in glob.glob(f'{file_path}/priors/ChIP_annotated/{tf}*annotation.txt'):
                chip_anno=pd.read_table(f)
                genel=chip_anno['Gene Name'].dropna().tolist()
                gene_count={k:genel.count(k) for k in genel}
        
                if tf not in prior_df.columns:
                    prior_df[f'{tf}']=[0]*prior_df.shape[0]
                    
                for k in gene_count:
                    if k not in prior_df.index:
                        prior_df.append({f'{tf}':1},ignore_index=True).fillna(0)
                    elif prior_df.loc[k,f'{tf}']==0:
                        prior_df.loc[k,f'{tf}']=1
                        
        prior_df.to_csv(f'{file_path}/priors/prior_atac_Miraldi_b_ChIP.tsv',sep='\t',header=True)


