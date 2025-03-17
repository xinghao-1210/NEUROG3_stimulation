library(glue)

# core_analysis_create_core_GRN.R
# Builds the following core GRNs from results of core TF GSEA:
# 1. core TFs + all its targets
# 2. core TFs + only its targets that are signature genes
# 3. core TFs + only its targets that are signature AND in a gene set of interest
print('--------------------------------------')
rm(list=ls())
options(stringsAsFactors=FALSE)

setwd('/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/HIOE_NEUROG3_induction')
#================== INPUTS ===================

# GRN file
maxcomb.path <- 'outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb'
maxcomb.par <- basename(maxcomb.path)
file_grn <- paste0(maxcomb.path,'/prior_atac_Miraldi_q_ChIP_bias10_maxComb_sp.tsv')
time.cond.list <- list(24,48,72,96) 

# FDR cutoff
fdr <- 0.05

for (time.cond in time.cond.list){
  # output directory
  dir_out <- glue('{maxcomb.path}/{time.cond}hpi_Cores')
  
  # GSEA file up
  file_gsea_up <- glue('{maxcomb.path}/GSEA/{maxcomb.par}_cut01_{time.cond}hpiSet_Praw0p1_dir_wCut0p0_minSet5/{time.cond}hpiSet_fdr100_up_adjp.txt')
  # GSEA file down
  file_gsea_down <- glue('{maxcomb.path}/GSEA/{maxcomb.par}_cut01_{time.cond}hpiSet_Praw0p1_dir_wCut0p0_minSet5/{time.cond}hpiSet_fdr100_down_adjp.txt')
  
  # Gene set enrichment file
  file_gene_set <- glue('inputs/geneSets/{time.cond}hpiSet.txt')
  
  # Gene subset of interest (e.g., cytokines, chemokines, receptors)
  # 'None' if not using
  file_subset <- 'None'
  # Gene set name for save filename
  id_subset <- 'None'
  
  # Gene set up & down names
  name_gene_set <- c(glue('{time.cond}hpi_up'), glue('{time.cond}hpi_down'))
  # Gene set name for save filename
  id_gene_set <- 'HIOE_NEUROG3'
  
  #=============================================
  
  # File save basename
  fileSave <- tools::file_path_sans_ext(substr(basename(file_grn), 1, nchar(basename(file_grn))-7))
  fileSave <- glue('Core_{fileSave}_fdr{fdr*100}_{id_gene_set}_{time.cond}hpiSet')
  
  # Load GRN
  GRN <- read.delim(file_grn, sep='\t', header=TRUE)
  
  # Load GSEA files
  # Up
  GSEAup <- read.delim(file_gsea_up, header=TRUE, sep='\t', row.names=1)
  # Down
  GSEAdown <- read.delim(file_gsea_down, header=TRUE, sep='\t', row.names=1)
  
  # Load gene set
  geneSet <- read.delim(file_gene_set, header=FALSE, row.names=2, sep=' ')
  
  # Over-expressed genes
  geneUp <- strsplit(geneSet[name_gene_set[1],2], split='[:|:]')
  geneUp <- as.character(geneUp[[1]])
  # Under-expressed genes
  geneDown <- strsplit(geneSet[name_gene_set[2],2], split='[:|:]')
  geneDown <- as.character(geneDown[[1]])
  geneAllDE <- union(geneUp,geneDown)
  
  uniqueTF <- unique(GRN$TF)
  
  # TFs whose targets are enriched in DE gene sets
  # Up
  subsetUpTF <- intersect(uniqueTF,colnames(GSEAup))
  idxKeepUp <- which(GSEAup[name_gene_set[1],subsetUpTF] < fdr)
  allUpTF <- subsetUpTF[idxKeepUp]
  # Down
  subsetDownTF <- intersect(uniqueTF,colnames(GSEAdown))
  idxKeepDown <- which(GSEAdown[name_gene_set[2],subsetDownTF] < fdr)
  allDownTF <- subsetDownTF[idxKeepDown]
  
  
  # 1. GRN: TFs with enriched targets in DE genes + all its targets
  allEnrichedTF <- union(allUpTF, allDownTF)
  idxKeepIntAll <- which(GRN$TF %in% allEnrichedTF)
  coreGRN_All <- GRN[idxKeepIntAll,]
  file_out <- file.path(dir_out, paste0(fileSave,'_All_sp.tsv'))
  write.table(coreGRN_All, file_out, sep='\t', quote=FALSE, row.names=FALSE)
  file_out <- file.path(dir_out, paste0(fileSave,'_All_gene_TF.txt'))
  writeLines(sort(unique(coreGRN_All$TF)), file_out)
  file_out <- file.path(dir_out, paste0(fileSave,'_All_gene_Target.txt'))
  writeLines(sort(unique(coreGRN_All$Target)), file_out)
  file_out <- file.path(dir_out, paste0(fileSave,'_All_gene_All.txt'))
  writeLines(sort(unique(union(coreGRN_All$TF,coreGRN_All$Target))), file_out)
  
  # 2. GRN: TFs with enriched targets in DE genes + only its targets that are DE genes
  idxKeepIntDE <- which(coreGRN_All$Target %in% geneAllDE)
  coreGRN_DE <- coreGRN_All[idxKeepIntDE,]
  file_out <- file.path(dir_out, paste0(fileSave,'_DE_sp.tsv'))
  write.table(coreGRN_DE, file_out, sep='\t', quote=FALSE, row.names=FALSE)
  file_out <- file.path(dir_out, paste0(fileSave,'_DE_gene_TF.txt'))
  writeLines(sort(unique(coreGRN_DE$TF)), file_out)
  file_out <- file.path(dir_out, paste0(fileSave,'_DE_gene_Target.txt'))
  writeLines(sort(unique(coreGRN_DE$Target)), file_out)
  file_out <- file.path(dir_out, paste0(fileSave,'_DE_gene_All.txt'))
  writeLines(sort(unique(union(coreGRN_DE$TF,coreGRN_DE$Target))), file_out)
  
  # 3. GRN: TFs with enriched targets in DE genes + only its targets that are DE genes and in gene set of interest
  if (file_subset != 'None'){
    # Load gene set
    subGene <- unique(readLines(file_subset))
    idxKeepSub <- which(coreGRN_DE$Target %in% subGene)
    coreGRN_DE_sub <- coreGRN_DE[idxKeepSub,]
    file_out <- file.path(dir_out, paste0(fileSave,'_DE_',id_subset,'_sp.tsv'))
    write.table(coreGRN_DE_sub, file_out, sep='\t', quote=FALSE, row.names=FALSE)	
    file_out <- file.path(dir_out, paste0(fileSave,'_DE_',id_subset,'_gene_TF.txt'))
    writeLines(sort(unique(coreGRN_DE_sub$TF)), file_out)
    file_out <- file.path(dir_out, paste0(fileSave,'_DE_',id_subset,'_gene_Target.txt'))
    writeLines(sort(unique(coreGRN_DE_sub$Target)), file_out)
    file_out <- file.path(dir_out, paste0(fileSave,'_DE_',id_subset,'_gene_All.txt'))
    writeLines(sort(unique(union(coreGRN_DE_sub$TF,coreGRN_DE_sub$Target))), file_out)
  }
}

#==============================================================================
print('--------------------------------------')
print('Done!')
