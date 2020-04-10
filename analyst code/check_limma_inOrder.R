# Check if differential expression results (limma output files)
# are sorted by adjusted p-value

# Read file
M3_results <- read.csv('3M_limma_results.csv',header=TRUE)
CLOT_results <- read.csv('CLOT_limma_results.csv', header=TRUE)
CHLR_results <- read.csv('CHLR_limma_results.csv', header=TRUE)

#output<-(dim(3M_results))


# Check if adjusted p value column is sorted
out_3M <- ifelse(!is.unsorted(M3_results$adj.P.Val),TRUE,FALSE)
out_CLOT <- ifelse(!is.unsorted(CLOT_results$adj.P.Val),TRUE,FALSE)
out_CHLR <- ifelse(!is.unsorted(CHLR_results$adj.P.Val),TRUE,FALSE)

output <- cbind(c('3M','CLOT','CHLR'),c(out_3M,out_CLOT,out_CHLR))

# write out the results to file
write.csv(output,'check_limmaSorted.csv')
