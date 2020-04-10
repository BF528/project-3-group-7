# Get top most significant 10 DE genes

# file name as command-line argument
args <- commandArgs(trailingOnly=TRUE)

# Read file
dat <- read.csv(args[1], header=TRUE)
dat <- dat[order(-dat$logFC),]

# Compute median
new_mat <- cbind(dat[2], dat$logFC)
#new_mat <- new_mat[order(new_mat$logFC),]


# write significant results to file
s_file<-paste0('top10DE_',args[1])
write.csv(new_mat,s_file)
