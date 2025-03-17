library("DESeq2") # Load DESeq2
library(ggplot2)
library(glue)

file = "/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/XZ_07-21-28-08-05-2020_RNA-seq_HIOE/Counts_TPM_mat/Genes_Counts_Matrix_rearranged.txt"        
# countData as imported count matrix with row-gene_ID column-sample_name (DESeq2 only accept int value)
database <- read.table(file, sep = "\t", header=T,row.names =1)
# using the names function to see names of the variables and which column of
# data to which they correspond
names(database)
# Select only comaprison group in Wells data
database <- database[,]
# set groups with condition (factor), the groups with same condition need to be together
condition <- factor(rep(c(rep("wt_0_24h",3),rep("wt_100_24h",3)),time=4))
# set "wt_0_24h" as reference level
condition <- relevel(condition, "wt_0_24h")
time <- factor(rep(c(24,48,72,96), each=6))
# colData (Dataframe) to assign condition to each group
coldata <- data.frame(row.names = colnames(database), condition, time)
# design as differential comparison matrix to indicate comparison gourps/conditions
# construct dds matrix based on countData, colData and design
dds <- DESeqDataSetFromMatrix(countData=database, colData=coldata, design=~condition)
head(dds)
# filter genes with at least one count
dds <- dds[ rowSums(counts(dds)) >= 1, ]

vsd <- vst(dds, blind=TRUE)
PCA_plot <- plotPCA(vsd, intgroup=c('time','condition'))
print(PCA_plot)

####Batch removal
# Batch removal
vsd <- vst(dds, blind=TRUE)
batch <- c(1:6)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch)
PCA_plot_limma <- plotPCA(vsd, intgroup=c('time', 'condition')) + ggtitle("PCA: HIOE")
print(PCA_plot_limma)

vsd.df <- as.data.frame(assay(vsd))
vsd.df <- cbind('Gene' = rownames(vsd.df), vsd.df)
rownames(vsd.df) <- NULL
write.table(vsd.df,quote=FALSE,
            file="/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/HIOE_NEUROG3_induction/inputs/geneExpression/RNAseq_24_DESeq2_VSDcounts.txt",sep='\t',row.names=F)

              