#Read in desq normalized count data for each MOA
ahr <- read.csv("ahr_deseq_norm_counts_ord.csv", header = TRUE)
car <- read.csv("car_deseq_norm_counts_ord.csv", header = TRUE)
cyto <- read.csv("cyto_deseq_norm_counts_ord.csv", header = TRUE)

#naming the column that has genes 
colnames(ahr)[1] <- "geneID"
colnames(car)[1] <- "geneID"
colnames(cyto)[2] <-"geneID"

#merging the genes that are common to the three MOA datasets into one dataframe 
merge_1 <-Reduce(merge, list(ahr, car, cyto[,2:5]))
merge_2 <-Reduce(merge, list(ahr[1:4], car[1:4], cyto[,2:5])) #excluding control treatments

#merge_1 <- merge(ahr[1:4], car[1:4], by = "geneID")
#merge_all <- merge(merge_1, cyto[,2:5], by= "geneID")

#making lables for the columns 
labels_1 <- c("ahr_7997", "ahr_7999", "ahr_8002", "ahr_control_1", "ahr_control_2", "ahr_control_3", "car_8020", "car_8036", "car_8046","car_control_1","car_control_2", "car_control_3","cyto_7987","cyto_7988","cyto_7989")

labels_2 <- c("ahr_7997", "ahr_7999", "ahr_8002", "car_8020", "car_8036", "car_8046", "cyto_7987","cyto_7988","cyto_7989")

#making heat map
merge_1_matrix <-as.matrix(merge_1[1:3000, 2:16])
View(heatmap(merge_1_matrix, labCol = labels_1))

#making heat map for the data frame without controls 
merge_2_matrix <-as.matrix(merge_2[1:3000, 2:10])
View(heatmap(merge_2_matrix, labCol = labels_2))

