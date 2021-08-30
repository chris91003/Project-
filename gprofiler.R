install.packages("gprofiler2")
library(gprofiler2)


#Load Df
df <- read.csv("kendall-liver-gene-clusters-high-variance.csv")
df

#subset column F from cluster

genes <- subset(df, select= c("F"))
genes


#subset Column A.2 from Cluser

A.2 <- subset(df, select= c("A.2"))
A.2


#convert column F to vector

gene_vector <- df$F
gene_vector
class(gene_vector)

class(genes)
View(df)

#convert column A.2 to vector

A.2_vector <- df$A.2
A.2_vector



#Gene Ontology Cluster F 

GO<- gconvert(query = c(gene_vector), organism = "hsapiens", 
              target="ENSG", mthreshold = Inf, filter_na = TRUE)

GO

View(GO)


#Gene Ontology Cluster A.2

A.2 <- gconvert(query = c(A.2_vector), organism = "hsapiens", 
              target="ENSG", mthreshold = Inf, filter_na = TRUE)
A.2



