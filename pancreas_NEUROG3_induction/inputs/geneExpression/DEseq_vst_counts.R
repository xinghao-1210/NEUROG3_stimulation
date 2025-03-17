library("DESeq2") # Load DESeq2
library(ggplot2)
library(glue)

file = "/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/XZ_03-04-08-12-2020_RNA-seq_Pancreas/Counts_TPM_mat/Genes_Counts_Matrix_rearranged.txt"        
# countData as imported count matrix with row-gene_ID column-sample_name (DESeq2 only accept int value)
database <- read.table(file = file, sep = "\t", header=T,row.names =1)
# using the names function to see names of the variables and which column of
# data to which they correspond
names(database)
# duplicate 0_8hpi as 0_0hpi data at 0 time point
database.0hpi<- database[,rep(1:3,2)]
database.merge <- merge(database.0hpi,database,by=0)
database <- database.merge[-1]
row.names(database) <- database.merge$Row.names
colnames(database)[1:9] <- c('wt_0_8h_0hpi_1','wt_0_8h_0hpi_2','wt_0_8h_0hpi_3','wt_100_8h_0hpi_1','wt_100_8h_0hpi_2','wt_100_8h_0hpi_3','wt_0_8h_8hpi_1','wt_0_8h_8hpi_2','wt_0_8h_8hpi_3')
# set groups with condition (factor), the groups with same condition need to be together
condition <- factor(rep(c(rep("wt_0_8h",3),rep("wt_100_8h",3)),time=5))
# set "wt_0_8h" as reference level
condition <- relevel(condition, "wt_0_8h")
time <- factor(rep(c(0,8,24,48,72), each=6))
# colData (Dataframe) to assign condition to each group
coldata <- data.frame(row.names = colnames(database),condition,time)
# design as differential comparison matrix to indicate comparison gourps/conditions
# construct dds matrix based on countData, colData and design
dds_tc <- DESeqDataSetFromMatrix(countData=database, colData=coldata, design=~condition+time+condition:time)
head(dds_tc)
# filter genes with at least one count
dds_tc <- dds_tc[ rowSums(counts(dds_tc)) >= 1, ]
dds_tc

vsd <- vst(dds_tc, blind=TRUE)
pcaData <- plotPCA(vsd, intgroup=c('time','condition'),returnData=T)
pcaData
plotPCA(vsd, intgroup=c('time','condition'))

####Batch removal
vsd <- vst(dds_tc, blind=TRUE)
batch <- c(1:3,1:3,1:3,1:3,5:7,5:7,8:10,8:10,5:7,5:7)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch)
plotPCA(vsd, intgroup=c('time','condition')) + ggtitle("PCA: pancreas")

vsd.df <- as.data.frame(assay(vsd))
vsd.df <- cbind('Gene' = rownames(vsd.df), vsd.df)
rownames(vsd.df) <- NULL
#vsd.df.100.8 <- vsd.df[ , c(1,11,12,13,17,18,19,23,24,25,29,30,31)]
write.table(vsd.df[-c(2:7) ],quote=FALSE,
            file="/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/pancreas_NEUROG3_induction/inputs/geneExpression/RNAseq_8_DESeq2_VSDcounts.txt",sep='\t',row.names=F)

              