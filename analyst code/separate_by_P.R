# Divide differential expression results (limma output file)
# into 2 files based on if p-adjust < 0.05 (significant).

# file name as command-line argument
args <- commandArgs(trailingOnly=TRUE)

# Read file
results <- read.csv(args[1], header=TRUE)
#results <- read.csv('example_limma_results.csv', header=TRUE)

# Subset of p-adjust that is significant
significant <- results[results$adj.P.Val < 0.05,]

# write significant results to file
s_file<-paste0('signif_',args[1])
write.csv(significant,s_file)

# Subset of p-adjust that is not significant
insignificant <- results[results$adj.P.Val >= 0.05,]

# Write insignificant results to file
is_file<-paste0('insignif_',args[1])
write.csv(insignificant,is_file)
