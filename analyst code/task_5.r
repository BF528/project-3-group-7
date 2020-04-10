# TASK 5
# Divide differential expression results (limma output file)
# into 2 files based on if p-adjust < 0.05 (significant).

# 3M
# Read file
signif <- read.csv('C:/Users/reva_/OneDrive/BU/Term 2/BF528/Projects/Project3/signif_3M_limma_results.csv', header=TRUE)

# Create jpeg file for graph output
#jpeg("task5_histogram.jpg")

# Create histogram of fold-change values from significant genes
hist(signif$logFC, main="3ME analysis",
     xlab="log(fold-change)")

# Create scatterplot of fold-change vs. nominal p-value
plot(signif$logFC, signif$P.Value, 
     main="3ME analysis",
     xlab="log(fold-change)", ylab="Nominal P-value")

# CLOT
# Read file
signif <- read.csv('C:/Users/reva_/OneDrive/BU/Term 2/BF528/Projects/Project3/signif_CLOT_limma_results.csv', header=TRUE)

# Create jpeg file for graph output
#jpeg("task5_histogram.jpg")

# Create histogram of fold-change values from significant genes
hist(signif$logFC, main="CLO analysis",
     xlab="log(fold-change)")

# Create scatterplot of fold-change vs. nominal p-value
plot(signif$logFC, signif$P.Value, 
     main="CLO analysis",
     xlab="log(fold-change)", ylab="Nominal P-value")

# CHLR
# Read file
signif <- read.csv('C:/Users/reva_/OneDrive/BU/Term 2/BF528/Projects/Project3/signif_CHLR_limma_results.csv', header=TRUE)

# Create jpeg file for graph output
#jpeg("task5_histogram.jpg")

# Create histogram of fold-change values from significant genes
hist(signif$logFC, main="CHL analysis",
     xlab="log(fold-change)")

# Create scatterplot of fold-change vs. nominal p-value
plot(signif$logFC, signif$P.Value, 
     main="CHL analysis",
     xlab="log(fold-change)", ylab="Nominal P-value")

# TASK 6
library("ggplot2")
# Manually load data
# Concordance scatter plot
concord <- as.data.frame(concord)
# make plot
ggplot(concord,aes(fill=position,x=chemical,y=value))+
  geom_bar(position="dodge",stat="identity")+
  scale_fill_brewer(palette="Reds")+
  labs(title="Concordance for Toxgroup1 DE genes",
       x ="Chemical", y = "Concordance")

# Scatter plot: Concordance vs number of DE genes
# load data
plots <- read.csv('C:/Users/reva_/OneDrive/BU/Term 2/BF528/Projects/Project3/plots.csv', header=TRUE)

# Microarray plot
plot(plots$ov_percent~plots$num_mic,
     main="Microarray", xlab="Number of sigificant genes",
     ylab="Concordance of DEGs", ylim=c(0,80))
text(x = plots$num_mic,
     y = plots$ov_percent,
     labels = plots$X, pos=3)

# RNA-seq plot
plot(plots$ov_percent~plots$num_rnaseq,
     main="RNA-seq", xlab="Number of sigificant genes",
     ylab="Concordance of DEGs", ylim=c(0,80))
text(x = plots$num_rnaseq,
     y = plots$ov_percent,
     labels = plots$X, pos=3)


