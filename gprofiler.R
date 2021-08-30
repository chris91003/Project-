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

F <- gconvert(query = c(gene_vector), organism = "hsapiens", 
              target="ENSG", mthreshold = Inf, filter_na = TRUE)

F

View(F)


#Gene Ontology Cluster A.2

A.2 <- gconvert(query = c(A.2_vector), organism = "hsapiens", 
              target="ENSG", mthreshold = Inf, filter_na = TRUE)
A.2



View(F)
View(A.2)


#Create dataframe that contains only gene names and function for each cluster
F_subset <- subset(F, select= c("name", "description"))
F_subset

A.2_subset <- subset(A.2, select= c("name", "description"))
A.2_subset


library(dplyr)


mergedDf <- F_subset %>% full_join(A.2_subset)

mergedDf

View(mergedDf)














#More Gene Ontology 

GO_A.2 <- gost(query = c(A.2), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL)

class(GO_A.2)

GO_A.2


GO_F <- gost(query = c(genes), 
                 organism = "hsapiens", ordered_query = FALSE, 
                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                 measure_underrepresentation = FALSE, evcodes = TRUE, 
                 user_threshold = 0.05, correction_method = "g_SCS", 
                 domain_scope = "annotated", custom_bg = NULL, 
                 numeric_ns = "", sources = NULL)



View(gostres2)

gostplot(genes, capped = TRUE, interactive= TRUE)

gostplot(A.2, capped= TRUE, interactive = TRUE )

