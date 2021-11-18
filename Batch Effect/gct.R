library(R.utils)
library(data.table)
library(ggplot2)
library(Hmisc)
library(dplyr)


#download.file("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
              destfile = "gtex_read_counts.gct.gz")


#gunzip("gtex_read_counts.gct.gz", remove= T)

data <- fread("gtex_read_counts.gct",data.table = FALSE)
View(data)




#Separate to skin and combine skin dataframe and read data frames

library(readr)


sample <- read_delim("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")


sample1 <- as.data.frame(sample)


skin1 <- sample1[which(sample1$SMTSD== 'Skin - Not Sun Exposed (Suprapubic)'), ]
View(skin1)

columns_rows <- colnames(data) %in% skin1$SAMPID
skin_data <- data[,columns_rows]

dat <- skin1[columns_rows,]

#Isolate gene Id from read data and combine to skin_data


Gene_ID <- data$Description
gene_id <- as.data.frame(Gene_ID)

View(gene_id)

final_skin <- cbind(gene_id, skin_data)
View(final_skin)

#Generate histogram of rows
hist(as.numeric(skin_data[2,]), xlab= "WASH7P", main= "WASH7P Histogram Read Count")
hist(as.numeric(skin_data[18,]), xlab= "RP11-34P13.18", main= "RP11-34P13.18 Histogram Read Count")

#Get mean of all rows

means <- rowMeans(skin_data)
means_column <- as.data.frame(means)
means_column



#append 

final_skin1 <- cbind(means_column, final_skin)
View(final_skin1)



#Histogram Visualization
hist(log10(means), breaks=100, xlab = "Read Counts", main= "")
density <- density(log10(means))


#Filters all data
nrow(final_skin)
skincounts <- as.numeric(final_skin)

filter <- rowSums(skin_data) > 650000


skincounts <- skin_data[filter,]

View(skincounts)


#Another way to filter with gene names
filter1 <- rowSums(skin_data[, 2:604]) > 650000
skincounts1 <- final_skin[filter1, ]
View(skincounts1)
class(skincounts1)

write.csv(skincounts1,"C:\\Users\\cdima\\OneDrive\\Desktop\\Project\\Batch Effect\\skincounts.csv", row.names = FALSE)

#Correlation

cor2 <- cor(t(skincounts1[, 2:605]), method= "spearman")
cor2 <- as.data.frame(cor2)
class(cor2)

#pastes row on top and column on end
skin_gene <- rbind(skincounts1$Gene_ID, cor2)
cor2$gene_id <- paste(skincounts1$Gene_ID)

#Pastes column on left
skin_gene1 <- cbind(skincounts1$Gene_ID, cor2)
skin_gene2 <- rbind(skincounts1$Gene_ID, cor2)

View(cor2)

cor1 <- cor(t(skincounts), method="spearman")
as.data.frame(cor1)
View(cor1)

skin_gene3 <- cbind(skincounts1$Gene_ID, skin_gene2)


cor2 <- cor1
cor2 <- as.data.frame(cor2)

write.csv(cor1,"C:\\Users\\cdima\\OneDrive\\Desktop\\Project\\newcor.csv", row.names = FALSE)
write.csv(skin_gene1,"C:\\Users\\cdima\\OneDrive\\Desktop\\Project\\Batch Effect\\genecorr.csv", row.names = FALSE)


#plot correlations Smad2 and smad4
plot(cor2$`46675`, cor2$`46742`, xlab = "SMAD2", ylab= "SMAD4")
cor(cor2$`46675`, cor2$`46742`)
summary(lm(cor2$`46675`~cor2$`46742`))
abline(lm(cor2$`46675`~cor2$`46742`))

#Merge the correlation data and geneid

cor2$Gene_ID <- NA
merge <- merge(gene_id, cor2, all.y=TRUE)


merge1 <- cor2[gene_id$Gene_ID,]

#12687 POLR2B TSHB 2490
#BIRC5	AURKB

cor(final_skin$Gene_ID== "POLR2B", final_skin$Gene_ID== "TSHB")
cor(final_skin$Gene_ID== "BIRC5", final_skin$Gene_ID == "AURKB")

row12687 <- skincounts[12687, ]
row2490 <- skincounts[2490, 1]

#Histogram filtered data
filtered_means <- rowMeans(skincounts)
hist(log10(filtered_means), breaks=100, xlab = "Filtered Read Counts", main= "")
filtered_density <- density(log10(filtered_means))


#Plot Density graphs
plot(density, col=1)
lines(filtered_density, col=2)


#Regression by gene

plot(row12687, row2490)

#Load in pos and neg corr genes
library(readr)
pos_corr <- read_delim("skin_pos_corr_genes.txt", delim = "\t", escape_double = FALSE, 
                                  col_names = FALSE, trim_ws = TRUE)

pos_corr <- as.data.frame(pos_corr)

new_pos <- pos_corr[pos_corr$X5 > 0.90,]
View(new_pos)

#Correlation Matrix and Corrleation matrix of pvalues
cor <- cor(skincounts)
spearman <- cor(skincounts, method = "pearson")

skincount <- as.matrix(skincounts)

pval <- rcorr(skincount, type = c("spearman"))

cor <-as.matrix(cor)
cor <- as.data.frame(cor)
class(cor)

write.csv(cor,"C:\\Users\\cdima\\OneDrive\\Desktop\\Project\\cor.csv", row.names = FALSE)
write.csv(spearman,"C:\\Users\\cdima\\OneDrive\\Desktop\\Project\\spearman.csv", row.names = FALSE)

#Correlation matrix of raw readss
library(corrplot)
source("http://www.sthda.com/upload/rquery_cormat.r")
cormat<-rquery.cormat(skincounts, graphType="heatmap")




skin2 <- skin_data

for(i in 1:nrow(skin_data)) {
  hist(as.numeric(skin_data[i,]))
}


for(i in 1:5){
  print(i)
}

#Remove unused samples
filter_data <- dplyr::filter(.data = skin1, SAMPID %in% colnames(data))
View(filter_data)

filt <- filter_data[columns_rows,]


#Normalization
normal <- DESeqDataSetFromMatrix(countData = skin_data,
                       colData = filter_data,
                       design = ~ 1)
results <- DESeq(normal)



#Volcano Plot
with(results, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))




#Load phenotype data and merge with Ischemic time
phenotype <- read.csv("gtex_phenotypes_v8.csv")

GTEx_SAMPID_to_SUBJID <- function(sampids){
  stringr::str_extract(string = sampids, pattern = "GTEX-[[:alnum:]]*")
}

phenotype_SUBJID <- GTEx_SAMPID_to_SUBJID(phenotype$SUBJID)
phenotype_SUBJID

skincounts2 <- skincounts
new <- GTEx_SAMPID_to_SUBJID(colnames(skincounts2))
colnames(skincounts2) <- new
View(skincounts2)

sampleid <-GTEx_SAMPID_to_SUBJID(colnames(skincounts))
class(sampleid)


col <- colnames(skincounts2) %in% phenotype$SUBJID
new_skincounts <- skincounts2[,col]

join <- left_join(new_skincounts, phenotype, by= "SUBJID")


#Remove first column and then get all samples from read count data into phenotype data

phenotype <- phenotype[,-1 ]
row <- colnames(skincounts2) %in% phenotype$SUBJID

row1 <- phenotype$SUBJID %in% colnames(skincounts2)


df2 <- phenotype[row1, ]


#Now have a dataframe saved as df2 which has the skin count data phenotypes

ischemic <- df2$TRISCHD
pca <- princomp(df2$TRISCHD, cor= TRUE, score= TRUE)       
summary(pca)
biplot(pca)
