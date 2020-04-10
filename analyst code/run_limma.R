library(limma)

# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_1_mic_info.csv',as.is=TRUE, header=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
  sep='\t',
  as.is=TRUE,
  header=TRUE,
  row.names=1,
)


## 3M
# subset the full expression matrix to just those in this comparison
rma.subset <- rma[paste0('X',samples$array_id[samples$chemical=='3-METHYLCHOLANTHRENE'|samples$chemical=='Control'])]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control','3-METHYLCHOLANTHRENE')
  )
)
colnames(design) <- c('Intercept','3-METHYLCHOLANTHRENE')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(t,'3M_limma_results.csv')


## CLOT
# subset the full expression matrix to just those in this comparison
rma.subset <- rma[paste0('X',samples$array_id[samples$chemical=='CLOTRIMAZOLE'|samples$chemical=='Control'])]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control','CLOTRIMAZOLE')
  )
)
colnames(design) <- c('Intercept','CLOTRIMAZOLE')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(t,'CLOT_limma_results.csv')


## CHLR
# subset the full expression matrix to just those in this comparison
rma.subset <- rma[paste0('X',samples$array_id[samples$chemical=='CHLOROFORM'|samples$chemical=='Control'])]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control','CHLOROFORM')
  )
)
colnames(design) <- c('Intercept','CHLOROFORM')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(t,'CHLR_limma_results.csv')
