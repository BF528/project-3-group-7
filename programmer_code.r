library(ggplot2)
library(reshape2)
library(DESeq2)
#set up matrix of counts
SRR1177987 = read.table("counts/SRR1177987_Aligned", sep = "")
SRR1177987 = SRR1177987[-1,]
SRR1177987 = SRR1177987[, 7]
SRR1177987 <- as.numeric(levels(SRR1177987))[SRR1177987]

SRR1177988 = read.table("counts/SRR1177988_Aligned", sep = "")
SRR1177988 = SRR1177988[-1,]
SRR1177988 = SRR1177988[, 7]
SRR1177988 <- as.numeric(levels(SRR1177988))[SRR1177988]

SRR1177989 = read.table("counts/SRR1177989_Aligned", sep = "")
SRR1177989 = SRR1177989[-1,]
SRR1177989 = SRR1177989[, 7]
SRR1177989 <- as.numeric(levels(SRR1177989))[SRR1177989]

SRR1177997 = read.table("counts/SRR1177997_Aligned", sep = "")
SRR1177997 = SRR1177997[-1,]
SRR1177997 = SRR1177997[, 7]
SRR1177997 <- as.numeric(levels(SRR1177997))[SRR1177997]

SRR1178002 = read.table("counts/SRR1178002_Aligned", sep = "")
SRR1178002 = SRR1178002[-1,]
SRR1178002 = SRR1178002[, 7]
SRR1178002 <- as.numeric(levels(SRR1178002))[SRR1178002]

SRR1178020 = read.table("counts/SRR1178020_Aligned", sep = "")
SRR1178020 = SRR1178020[-1,]
SRR1178020 = SRR1178020[, 7]
SRR1178020 <- as.numeric(levels(SRR1178020))[SRR1178020]

SRR1178036 = read.table("counts/SRR1178036_Aligned", sep = "")
SRR1178036 = SRR1178036[-1,]
SRR1178036 = SRR1178036[, 7]
SRR1178036 <- as.numeric(levels(SRR1178036))[SRR1178036]

SRR1178046 = read.table("counts/SRR1178046_Aligned", sep = "")
SRR1178046 = SRR1178046[-1,]
SRR1178046 = SRR1178046[, 7]
SRR1178046 <- as.numeric(levels(SRR1178046))[SRR1178046]

SRR1177999 = read.table("counts/SRR1177999_Aligned", sep = "")
SRR1177999 = SRR1177999[-1,]
SRR1177999 = SRR1177999[, 7]
SRR1177999 <- as.numeric(levels(SRR1177999))[SRR1177999]

geneid = read.table("counts/SRR1177987_Aligned", sep = "")
geneid = geneid[, 1]
geneid = geneid[-1]


counts = data.frame(geneid, SRR1177987, SRR1177988, SRR1177989, SRR1177997, SRR1177999, SRR1178002, SRR1178020, SRR1178036, SRR1178046)
write.csv(counts, "counts.csv")

colnames(counts)
row.names(counts) = counts[, 1]
counts = counts[,-1]


#boxplot for counts
ggplot(melt(counts), aes(variable, log(value))) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Samples", y="log(Counts)") +
  theme_classic() +
  ylim(0, 12) +
  theme(axis.text.x = element_text(angle=30, size=10, hjust=1), 
        axis.title.x = element_text(size = 13, face="bold"),
        axis.title.y = element_text(size = 13, face="bold")) 

#load in control samples
control_counts = read_csv("control_counts.csv")
control_counts = as.data.frame(control_counts)
row.names(control_counts) = control_counts[,1]
control_counts = control_counts[,-1]

all_counts = cbind(counts, control_counts)

#splitting up groups
mdata = read.csv('toxgroups/toxgroup_1_rna_info.csv')

ahr_cmc = mdata[mdata$mode_of_action=="AhR" | mdata$vehicle=="CMC_.5_%",]
car = rbind(mdata[mdata$mode_of_action=="CAR/PXR",], mdata[mdata$mode_of_action=="Control" & mdata$vehicle=="CORN_OIL_100_%",])
cyto = rbind(mdata[mdata$mode_of_action=="Cytotoxic",], mdata[mdata$mode_of_action=="Control" & mdata$vehicle=="CORN_OIL_100_%",])

ahr_cnts = cbind(all_counts$SRR1177997, all_counts$SRR1177999, all_counts$SRR1178002, all_counts$SRR1178000, all_counts$SRR1178005, all_counts$SRR1178007)
row.names(ahr_cnts) <- rownames(all_counts)
ahr_cnts = as.data.frame(ahr_cnts)
colnames(ahr_cnts) <- ahr_cmc$Run

car_cnts = cbind(all_counts$SRR1178020, all_counts$SRR1178036, all_counts$SRR1178046, all_counts$SRR1177973, all_counts$SRR1178016, all_counts$SRR1178019)
row.names(car_cnts) <- rownames(all_counts)
car_cnts = as.data.frame(car_cnts)
colnames(car_cnts) <- car$Run

cyto_cnts = cbind(all_counts$SRR1177987, all_counts$SRR1177988, all_counts$SRR1177989, all_counts$SRR1177973, all_counts$SRR1178016, all_counts$SRR1178019)
row.names(cyto_cnts) <- rownames(all_counts)
cyto_cnts = as.data.frame(cyto_cnts)
colnames(cyto_cnts) <- cyto$Run

# filter out rows that have any zeros for funzies
ahr_cnts <- subset(ahr_cnts,rowSums(ahr_cnts==0)==0)
car_cnts <- subset(car_cnts,rowSums(car_cnts==0)==0)
cyto_cnts <- subset(cyto_cnts,rowSums(cyto_cnts==0)==0)

##differential expression (repeated for each group)
##AHR
#create the DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = ahr_cnts,
  colData = ahr_cmc,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')

# run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c('mode_of_action','AhR','Control'))
res <- lfcShrink(dds, coef=2)

# write out DE results
write.csv(res,'ahr_deseq_results.csv')

# write out matrix of normalized counts
write.csv(counts(dds,normalized=TRUE),'ahr_deseq_norm_counts.csv')
##CAR
#create the DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = car_cnts,
  colData = car,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')

# run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c('mode_of_action','CAR/PXR','Control'))
res <- lfcShrink(dds, coef=2)

# write out DE results
write.csv(res,'car_deseq_results.csv')

# write out matrix of normalized counts
write.csv(counts(dds,normalized=TRUE),'car_deseq_norm_counts.csv')
##CYTO
#create the DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = cyto_cnts,
  colData = cyto,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')

# run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c('mode_of_action','Cytotoxic','Control'))
res <- lfcShrink(dds, coef=2)

# write out DE results
write.csv(res,'cyto_deseq_results.csv')

# write out matrix of normalized counts
write.csv(counts(dds,normalized=TRUE),'cyto_deseq_norm_counts.csv')

#reload results
ahr_deseq <- read.csv("DESeq/ahr_deseq_results.csv")
ahr_deseq <- as.data.frame(ahr_deseq)
row.names(ahr_deseq) <- ahr_deseq[, 1]
ahr_deseq <- ahr_deseq[,-1]

car_deseq <- read.csv("DESeq/car_deseq_results.csv")
car_deseq <- as.data.frame(car_deseq)
row.names(car_deseq) <- car_deseq[, 1]
car_deseq <- car_deseq[,-1]

cyto_deseq <- read.csv("DESeq/cyto_deseq_results.csv")
cyto_deseq <- as.data.frame(cyto_deseq)
row.names(cyto_deseq) <- cyto_deseq[, 1]
cyto_deseq <- cyto_deseq[,-1]

#sort p-values
ahr_deseq_padj <- ahr_deseq[order(ahr_deseq$padj),]
ahr_deseq_padj <- ahr_deseq_padj[ahr_deseq_padj$padj < 0.05,]
ahr_deseq_padj <- ahr_deseq_padj[complete.cases(ahr_deseq_padj), ]

car_deseq_padj <- car_deseq[order(car_deseq$padj),]
car_deseq_padj <- car_deseq_padj[car_deseq_padj$padj < 0.05,]
car_deseq_padj <- car_deseq_padj[complete.cases(car_deseq_padj), ]

cyto_deseq_padj <- cyto_deseq[order(cyto_deseq$padj),]
cyto_deseq_padj <- cyto_deseq_padj[cyto_deseq_padj$padj < 0.05,]
cyto_deseq_padj <- cyto_deseq_padj[complete.cases(cyto_deseq_padj), ]

#create padj filtered csv files
write.csv(ahr_deseq_padj, file="DESeq/ahr_deseq_results_padj.csv")
write.csv(car_deseq_padj, file="DESeq/car_deseq_results_padj.csv")
write.csv(cyto_deseq_padj, file="DESeq/cyto_deseq_results_padj.csv")

#histograms
ggplot(ahr_deseq_padj, aes(ahr_deseq_padj$log2FoldChange)) +
  geom_histogram(bins = 50, binwidth = 0.5, color="black") +
  xlab("log2 Fold Change") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 13, face="bold"),
        axis.title.y = element_text(size = 13, face="bold"))
  
ggplot(car_deseq_padj, aes(car_deseq_padj$log2FoldChange)) +
  geom_histogram(bins = 50, binwidth = 0.5, color="black") +
  xlab("log2 Fold Change") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 13, face="bold"),
        axis.title.y = element_text(size = 13, face="bold"))

ggplot(cyto_deseq_padj, aes(cyto_deseq_padj$log2FoldChange)) +
  geom_histogram(bins = 50, binwidth = 0.5, color="black") +
  xlab("log2 Fold Change") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 13, face="bold"),
        axis.title.y = element_text(size = 13, face="bold"))

#volcano plots
av = ahr_deseq[order(ahr_deseq$log2FoldChange),]
cv = car_deseq[order(car_deseq$log2FoldChange),]
cy = cyto_deseq[order(cyto_deseq$log2FoldChange),]

ggplot(ahr_deseq, aes(log2FoldChange, -log(pvalue))) + 
  geom_point() +
  geom_hline(yintercept = -log(0.05), colour="red", linetype="dashed") +
  geom_point(data=tail(av, 10), col="red") +
  geom_point(data=head(av, 10), col="dodgerblue2") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 13, face="bold"),
        axis.title.y = element_text(size = 13, face="bold"))
  
ggplot(car_deseq, aes(log2FoldChange, -log(pvalue))) + 
  geom_point() +
  geom_hline(yintercept = -log(0.05), colour="red", linetype="dashed") +
  geom_point(data=tail(cv, 10), col="red") +
  geom_point(data=head(cv, 10), col="dodgerblue2") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 13, face="bold"),
        axis.title.y = element_text(size = 13, face="bold"))

ggplot(cyto_deseq, aes(log2FoldChange, -log(pvalue))) +
  geom_point() +
  geom_hline(yintercept = -log(0.05), colour="red", linetype="dashed") +
  geom_point(data=tail(cy, 10), col="red") +
  geom_point(data=head(cy, 10), col="dodgerblue2") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 13, face="bold"),
        axis.title.y = element_text(size = 13, face="bold"))

#reload normalized counts
ahr_deseq_cnts <- read.csv("DESeq/normalized_counts/ahr_deseq_norm_counts.csv")
ahr_deseq_cnts <- as.data.frame(ahr_deseq_cnts)
row.names(ahr_deseq_cnts) <- ahr_deseq_cnts[, 1]
ahr_deseq_cnts <- ahr_deseq_cnts[,-1]

car_deseq_cnts <- read.csv("DESeq/normalized_counts/car_deseq_norm_counts.csv")
car_deseq_cnts <- as.data.frame(car_deseq_cnts)
row.names(car_deseq_cnts) <- car_deseq_cnts[, 1]
car_deseq_cnts <- car_deseq_cnts[,-1]

cyto_deseq_cnts <- read.csv("DESeq/normalized_counts/cyto_deseq_norm_counts.csv")
cyto_deseq_cnts <- as.data.frame(cyto_deseq_cnts)
row.names(cyto_deseq_cnts) <- cyto_deseq_cnts[, 1]
cyto_deseq <- cyto_deseq_cnts[,-1]

ahr_deseq_cnts <- ahr_deseq_cnts[match(row.names(ahr_deseq), row.names(ahr_deseq_cnts)),]
car_deseq_cnts <- car_deseq_cnts[match(row.names(car_deseq), row.names(car_deseq_cnts)),]
cyto_deseq_cnts <- cyto_deseq_cnts[match(row.names(cyto_deseq), row.names(cyto_deseq_cnts)), ]

write.csv(ahr_deseq_cnts, "DESeq/normalized_counts/ahr_deseq_norm_counts_ord.csv")
write.csv(car_deseq_cnts, "DESeq/normalized_counts/car_deseq_norm_counts_ord.csv")
write.csv(cyto_deseq_cnts, "DESeq/normalized_counts/cyto_deseq_norm_counts_ord.csv")