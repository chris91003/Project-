install.packages("gprofiler2")
library(gprofiler2)


#Load Df
df <- read.csv("kendall-liver-gene-clusters-high-variance.csv")
df

#subset columns from cluster

genes <- subset(df, select= c("F"))
genes


#convert column to vector

gene_vector <- df$F
gene_vector
class(gene_vector)

class(genes)
View(df)



#Gene Ontology 

GO<- gconvert(query = c(gene_vector), organism = "hsapiens", 
              target="ENSG", mthreshold = Inf, filter_na = TRUE)

GO

View(GO)
