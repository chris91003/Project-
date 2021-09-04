#Iterate through data frame and take each column and convert to Vector 
#iterate through columns and maybe see if there are GO terms that are in multiple clusters?
#iterate using other tissues????
library(gprofiler2)

#load in data
df <- read.csv("kendall-liver-gene-clusters-high-variance.csv")
df

View(df)


#for loop to iterate through df and extract columns and convert to vector 

for(i in 1:length(df)) {
  column <- c(df, i)
}



#Now have list of vectors that can be inserted into Gprofiler
cluster_A.2 <- gconvert(query = c(column$X), organism = "hsapiens", 
                        target="ENSG", mthreshold = Inf, filter_na = TRUE)
cluster_A.2



