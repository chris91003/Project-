#zscore to represent how much a sample is associated with a cluster
# if z score of 3 that means genes had very high expression relative to other
#each sample is tissue 



install.packages("tidyverse")
library("tidyverse")

install.packages("dplyr")
library("dplyr")

library(Hmisc)
library(ggplot2)



getwd()
df <- read.csv("kendall-liver-cluster-profiles.csv")
phenotype <- read.csv("gtex_phenotypes_v8.csv")



View(df)
describe(df)



#Exploratory data an analysis to learn more about the data




#for loop to make histogram of all columms and see that it is left skewed and majority of data has higher gene expression

hist=for(col in 2:ncol(df)) {
  hist(df[,col],xlab="Distribution of Zscores", main= "Histogram for Zscores", 
       col='Blue')
}



# See if correlation is significant and Graph Correlations 
install.packages("Hmisc")
library("Hmisc")

corr_matrix <- rcorr(as.matrix(new_df))
corr_matrix

#Remove dataframe column and create a heatmap
new_df <- df[2:12]
new_df

corr= cor(new_df)

num_df <- data.matrix(new_df)
heatmap(num_df)



install.packages("corrplot")
library("corrplot")

corrplot(corr, method= "circle", addCoef.col = "black", title= "Correlation Matrix of Data",
         mar=c(0,0,2,0))



#Scatter plot and correlation of A.1 and A.2 column (looks like neg corr)


x <- new_df$A.1
y <- new_df$A.2
plot(x,y, xlab= "A.1 Zscore", ylab= "A.2 Zscore")
abline(lm(y~x), col="red")
cor(new_df$A.1,new_df$A.2)


#Boxplot of A1 and A2

par(mfrow=c(1,2))
boxplot(x, xlab= "Boxplot of A.1", ylab= "ZScore",col= "red")
boxplot(y, xlab= "Boxplot of A.2", ylab= "ZScore",col= "blue")


#Bar plot of A1 vs A2
bar <- ggplot(data= new_df, aes(x= A.1, y= A.2, fill= x, y))+ 
  geom_bar(stat= "identity", width= 0.01)+
  coord_flip()

bar    




#Finding highly correlated clusters
max_corr <- which(abs(corr) > 0.5 , arr.ind = TRUE)
max_corr


corr_df <- as.data.frame.matrix(corr)
corr_df
corr_max <- corr_df[corr_df$cor > 0.5 & corr_df$cor <0.9]
corr_max

which(corr_df > 0.5 & corr_df <0.9)

#Function to Sort into specific columns
index_function <- function(corr_df){
  index <- which(corr_df > 0.5 & corr_df <0.9, arr.ind = TRUE)
  return(index)}

index_function(corr_df)


#Finding negatively correlated clusters
negative_function <- function(corr_df){
  index <- which(corr_df < 0.01, arr.ind = TRUE)
  return(index)}

negative_function(corr_df)

#Now extract indices 17, 19, 47, 62, 65, and 82 and put into another dataframe
corr1 <- corr_df[corr_df > 0.5 & corr_df <0.9]
corr1


#Looks like specific clusters have high correlation, wonder why that is? 
#Maybe because these clusters have interactions with each other?



#Join zscore and phenotype data

GTEx_SAMPID_to_SUBJID <- function(sampids){
  ## A simple function to quickly turn SAMPIDs into SUBJIDs to connect
  ## individuals to their phenotype
  stringr::str_extract(string = sampids, pattern = "GTEX-[[:alnum:]]*")
}

zscore_SUBJID <- GTEx_SAMPID_to_SUBJID(df$X)
zscore_SUBJID

#Replace subject ID with new Subj ID and join
df$X <- zscore_SUBJID

#rename zscore data column
names(df)[1] <- "SUBJID"
join <- left_join(df, phenotype, by= "SUBJID")
View(join)

join$MHDRNKSTS <- ifelse(join$MHDRNKSTS == "Yes", 1, 0)



ggplot(join) + aes(MHDRNKSTS, A.1) +
  geom_jitter()
