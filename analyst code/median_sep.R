# Divide DE genes into above and below median files

# file name as command-line argument
args <- commandArgs(trailingOnly=TRUE)

# Read file
deseq <- read.csv(args[1], header=TRUE)

# Compute median
med <- median(deseq$baseMean)

# Subset of p-adjust that is significant
above <- deseq[deseq$baseMean>=med,]

# write significant results to file
s_file<-paste0('aboveMedian_',args[1])
write.csv(above,s_file)

# Subset of p-adjust that is not significant
below <- deseq[deseq$baseMean < med,]

# Write insignificant results to file
is_file<-paste0('belowMedian_',args[1])
write.csv(below,is_file)
