# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import re
import glob
import pandas as pd
import numpy as np

file_path='/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/pancreas_NEUROG3_induction/inputs'
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
#!cat *.bed | bedtools sort | bedtools merge > output.bed

#==============================================
# TFs to list, intersect with RNA-seq expression
# normGene: vsd normalized gene x sample matrix
# genesforTFA: ALL genes from RNA-seq across samples
# targetGene: Differentially expressed genes from RNA-seq across samples
# potRegs: Differentially expressed genes from RNA-seq across samples INTERSECT with motif2gene list
#==============================================
tfs_anno=pd.read_table(f'{file_path}/geneSets/human_TF_set.txt')
rna=pd.read_table(f'{file_path}/geneExpression/RNAseq_8_DESeq2_VSDcounts.txt')

tf_add=['NEUROG3','PDX1','NKX2-2','NKX6-1','FOXA1','FOXA2','MAFB','MAFA','RFX6','GATA4','GATA6','NEUROD1','NEUROD2','PAX4','PAX6']+[x.strip() for x in open(f'{file_path}/geneSets/Gehart2019_eec.txt')]
tf_add=set(tf_add)
targetGene=set()
# target gene as union of deg in both pancreas and intestine
for rna_exp in ['XZ_03-04-08-12-2020_RNA-seq_Pancreas']:
    for f in glob.glob(f'/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/{rna_exp}/Counts_TPM_mat/analysis/results/*hpi_filtered.txt'):
        de_pd=pd.read_table(f)
        de_gene=[x for x in de_pd.Gene.tolist() if x.isupper()]
        targetGene=targetGene | set(de_gene)
targetGene=(targetGene | tf_add) & set(rna.iloc[:,0].tolist())
with open(f'{file_path}/targRegLists/targetGenes_names.txt','w') as f:
    f.writelines('\n'.join(targetGene))
    
potRegs=list(set(tfs_anno.Gene.tolist()) & targetGene)  
with open(f'{file_path}/targRegLists/potRegs_names.txt','w') as f:
    f.writelines('\n'.join(potRegs))

# background (DEFAULT all genes, optionally slice relieved padj as background) and genesForTFA
genesForTFA=rna
genesForTFA.to_csv(f'{file_path}/targRegLists/genesForTFA.txt',sep='\t',header=True,index=False)
genesForTFA.iloc[:,0].to_csv(f'{file_path}/geneSets/background_set.txt',sep='\t',header=False,index=False)

# background with DE
for p in [0.05, 0.15, 0.25]:
    bg_gene=set() 
    for f in glob.glob('/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/XZ_03-04-08-12-2020_RNA-seq_Pancreas/Counts_TPM_mat/analysis/results/*hpi.txt'):
        de_pd=pd.read_table(f)
        de_gene=de_pd[(abs(de_pd.log2FoldChange)>1) & (de_pd.padj<p)].Gene.tolist()
        bg_gene=bg_gene | set(de_gene)
    with open(f'{file_path}/geneSets/background_set_de{p}.txt','w') as f:
        f.writelines('\n'.join(bg_gene))
#==============================================
# add annotated ChIP onto merged ATAC priors
#==============================================
tfl=['NEUROG3','PDX1','NKX2-2','NKX6-1','FOXA1','FOXA2','MAFB']
weight=1

prior_df_q=pd.read_csv(f'{file_path}/priors/prior_atac_Miraldi_q.tsv',sep='\t',index_col=0)
prior_df_b=pd.read_csv(f'{file_path}/priors/prior_atac_Miraldi_b.tsv',sep='\t',index_col=0)

for p, prior in enumerate([prior_df_q, prior_df_b]):
    prior_df=prior_df_q
    
    if p==0:
        for tf in tfl:
            for f in glob.glob(f'{file_path}/priors/ChIP_annotated/{tf}*annota*.txt'):
                chip_anno=pd.read_table(f)
                genel=chip_anno.groupby('Gene Name')['Peak Score'].count()
                gene_count=genel.to_dict()
                # genel=chip_anno['Gene Name'].dropna().tolist()
                # gene_count={k:genel.count(k) for k in genel}
    
                # only use ChIP-seq data for available TFs
                prior_df[f'{tf}']=[0]*prior_df.shape[0]
                
                for k in gene_count:
                    if k not in prior_df.index:
                        prior_df.loc[k]=[0]*prior_df.shape[1]
                        prior_df.loc[k,f'{tf}']=gene_count[k]*weight
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
        
                prior_df[f'{tf}']=[0]*prior_df.shape[0]
                    
                for k in gene_count:
                    if k not in prior_df.index:
                        prior_df.loc[k]=[0]*prior_df.shape[1]
                        prior_df.loc[k,f'{tf}']=1
                    elif prior_df.loc[k,f'{tf}']==0:
                        prior_df.loc[k,f'{tf}']=1
                        
        prior_df.to_csv(f'{file_path}/priors/prior_atac_Miraldi_b_ChIP.tsv',sep='\t',header=True)


# prior_df_q=pd.read_csv(f'{file_path}/priors/prior_atac_Miraldi_q_ChIP.tsv',sep='\t',index_col=0)
# prior_df_b=pd.read_csv(f'{file_path}/priors/prior_atac_Miraldi_b_ChIP.tsv',sep='\t',index_col=0)
# prior_df_q.loc['GIP','PDX1']


