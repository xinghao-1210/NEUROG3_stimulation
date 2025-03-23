# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import re
import pandas as pd
import numpy as np

#==============================================
# motifs division and motif_2_gene
#==============================================
with open('JASPAR2018_combined_matrices_pfm.txt','r') as f:
    motif=f.readlines()

motif_2_gene={}

for s in motif:
    if '>' in s:
        motif_sub=motif[motif.index(s)+1:motif.index(s)+5]
        motif_title=motif[motif.index(s)].split(' ')[0][1:]
        motif_gene=motif[motif.index(s)].split(' ')[1][:-1]
        if '(' and ')' in motif_gene:
            motif_gene=motif_gene.split('(')[0]
        motif_2_gene[f'{motif_title}']=motif_gene

        A=[]
        for l in motif_sub:
            nl=[float(x) for x in re.split(' |\t',l) if '.' in x]
            A.append(nl)
        A=np.asarray(A)
        A=A/np.sum(A,axis=0)
        np.savetxt(f'./motifs/{motif_title}.pfm', A, delimiter='\t')

motif_2_gene_pd=pd.DataFrame.from_dict(motif_2_gene, orient='index')
motif_2_gene_pd.to_csv('tbl_motif_2_gene.tsv',sep='\t',header=False)


#==============================================
# gtf to bed file as gene body
#==============================================
!gunzip -c gencode.v34lift37.annotation.gtf.gz |  awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$14,$7}}' | tr -d '";' | head
!gunzip -c gencode.v34lift37.annotation.gtf.gz |  awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$14,$7}}' | tr -d '";' > gene_body.bed


#==============================================
# merge ATAC peaks
#==============================================
!for f in *.narrowpeaks_peaks.txt;do sort -k1,1 -k2,2n $f > ${f/.narrowpeaks_peaks.txt/.sorted.narrowpeaks_peaks.txt}; done
!bedtools intersect -a ATAC-seq--WT-Pan-100-8hpi--D9--7-30--E01744.sorted.narrowpeaks_peaks.txt -b ATAC-seq--WT-Pan-100-8hpi--D9--7-26--E01720.sorted.narrowpeaks_peaks.txt ATAC-seq--WT-Pan-100-8h--D12-3dpi--7-30--E01740.sorted.narrowpeaks_peaks.txt ATAC-seq--WT-Pan-100-8h--D12-3dpi--7-26--E01714.sorted.narrowpeaks_peaks.txt > ../ATAC_peaks.bed


#==============================================
# add annotated ChIP onto merged ATAC priors
#==============================================
tf='NEUROG3'
chip_anno=pd.read_table('/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/pancreas_NEUROG3_induction/priorConstruction/inputs/ChIP_annotated/09-17ngn3_30M_peaksfdr0.05_annotation.txt')
genel=chip_anno['Gene Name'].dropna().tolist()
gene_count={k:genel.count(k) for k in genel}

for prior in ['prior_atac_q','prior_atac_b']:
    prior_df=pd.read_csv(f'/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/pancreas_NEUROG3_induction/priorConstruction/outputs/{prior}.tsv',sep='\t',index_col=0)
    if tf not in prior_df.columns:
        prior_df[f'{tf}']=[0]*prior_df.shape[0]
    if '_q' in prior:
        for k in gene_count:
            if k not in prior_df.index:
                prior_df.append({f'{tf}':gene_count[k]},ignore_index=True).fillna(0)
            else:
                prior_df.loc[k,f'{tf}']=prior_df.loc[k,f'{tf}']+gene_count[k]
        prior_df.to_csv(f'/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/pancreas_NEUROG3_induction/priorConstruction/outputs/{prior}_ChIP.tsv',sep='\t',header=True)
                                                         
    if '_b' in prior:
        for k in gene_count:
            if k not in prior_df.index:
                prior_df.append({f'{tf}':1},ignore_index=True).fillna(0)
            elif prior_df.loc[k,f'{tf}']==0:
                prior_df.loc[k,f'{tf}']=1
        prior_df.to_csv(f'/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/pancreas_NEUROG3_induction/priorConstruction/outputs/{prior}_ChIP.tsv',sep='\t',header=True)
        
        