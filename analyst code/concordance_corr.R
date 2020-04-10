# Compute concordance
# ARGS: conconrdance.R <limma.csv> <deseq.csv>

# Input files
args=commandArgs(trailingOnly=TRUE)

# microarray data
limma <- read.csv(args[1],header=TRUE)
# RNA-seq data
deseq2 <- read.csv(args[2],header=TRUE)
# Mapping reference
map_ref <- read.csv('/project/bf528/project_3/refseq_affy_map.csv',header=TRUE)

# filter
deseq2 <- deseq2[deseq2$padj<0.05,]
deseq2 <- na.omit(deseq2)

limma_match <- map_ref$PROBEID %in% limma$X
deseq_match <- map_ref$REFSEQ %in% deseq2$X
index <- limma$logFC[match(map_ref$PROBEID,limma$X)]
index2 <- deseq2$log2FoldChange[match(map_ref$REFSEQ,deseq2$X)]
test <- match(map_ref$PROBEID,limma$X)

same_FCdir <- sign(index)==sign(index2)

intsect_two <- limma_match & deseq_match & same_FCdir
itsct_num <- sum(intsect_two)
total <- sum(limma_match | deseq_match)
'un_c'
2*itsct_num/total

n0 <- itsct_num
n1 <- dim(deseq2)[1]
n2 <- dim(limma)[1]
N <- 25225

new_intersection <- (n0*N-n1*n2)/(n0+N-n1-n2)
'c'
2*new_intersection/total
