#Iterate through data frame and take each column and convert to Vector 
#iterate through columns and maybe see if there are GO terms that are in multiple clusters?
#iterate using other tissues????
library(gprofiler2)
library(dplyr)


install.packages("rentrez")
library(rentrez)

library(pylr)
library(ggplot2)

#load in data
df <- read.csv("kendall-liver-gene-clusters-high-variance.csv")
df

View(df)

df$A.1

#for loop to iterate through df and extract columns and convert to vector 

for(i in 1:ncol(df)) {
  column <- c(df, i)
}



#Multiple query list for Gost Function

GO <- gost(query= list("Column A.1"= c(column$A.1),
                 "Column A.2"= c(column$A.2),
                 "Column B"= c(column$B),
                  "Column C"= c(column$C),
                  "Column D"= c(column$D),
                  "Column E"= c(column$E),
                  "Column F"= c(column$F),
                  "Column G"= c(column$G),
                  "Column H"= c(column$H),
                  "Column I.1"= c(column$I.1),
                  "Column I.2"=c(column$I.2)),
evcodes = TRUE)
                        

#Splitting Resulting Dataframe
result <- GO$result

result
split1 <- split(result, f=result$query)
split
View(split1)

#Iterate over list and extract individual dataframes

#Extract dataframes one by one



A.1 <- split1$`Column A.1`
A.1
View(A.1)

class(A)
class(column$A.1)
View(column$A.1)



A.2 <- split$`Column A.2`

View(A.2)

lapply(split, )



#Look at individual genes and where they are located 

A.1_list <- A.1[1,16]

class(A.1_list)
entrez_dbs()


ge <-entrez_search(db="gene", term="Homo Sapiens", id=c(A.1_list))
View(ge)

#Now have list of vectors that can be inserted into Gprofiler
GO_F <- gost(query = list("Column F"= c(column$F)), 
             organism = "hsapiens", ordered_query = FALSE, 
             multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
             measure_underrepresentation = FALSE, evcodes = TRUE, 
             user_threshold = 0.05, correction_method = "g_SCS", 
             domain_scope = "annotated", custom_bg = NULL, 
             numeric_ns = "", sources = NULL)


GO_F


#Gprofiler Result Dataframe
results <- GO_F$result
class(results)
View(results)


#Tabular format of GO IDs and Intersecting Genes
intersect<- data.frame("Go_ID"= c(results$term_id), 
                       "Intersecting Genes"= c(results$intersection),
                       "Description" = c(results$term_name))

intersect

View(intersect)



#Dataframe Exploration









#Maybe do GO for each column and analyze separately?


#Comparing Positively Correlated Columns > 0.5
A.2_F <- gost(query= list("Column F"= c(column$F),
                          "Column A.2"= c(column$A.2)),
              evcodes = TRUE)

A.2_H <- gost(query= list("Column A.2"= c(column$A.2),
                          "Column H"= c(column$H)),
              evcodes = TRUE)

D_F <- gost(query= list("Column D"= c(column$D),
                          "Column F"= c(column$F)),
              evcodes = TRUE)


#Results for F and A.2
A.2_Fresults <- A.2_F$result
A.2_Fresults
View(A.2_Fresults)

View(D_F$result)

#Iterate over column to find duplicate GO Terms
for(i in 1:length(A.2_Fresults$term_id)){
  if(A.2_Fresults$term_id[i]== i )
}


#Dataframe Exploration
barplot(A.2_Fresults$term_size)


#Tabular format of GO IDs and Intersecting Genes for A2 and F cluster
A.2_data<- A.2_Fresults[1:27, c("query", "term_id","p_value", 
                                "intersection","term_name")]
F_data<- A.2_Fresults[28:38, c("query", "term_id", "p_value",
                               "intersection", "term_name")]

View(A.2_data)

#Generate Manhattan Plots 

A.2_Fplot <- gostplot(A.2_F, capped= TRUE, interactive = TRUE)

A.2_Hplot <- gostplot(A.2_H, capped= TRUE, interactive = TRUE)

D_Fplot <- gostplot(D_F, capped= TRUE, interactive = TRUE)

A.2_Hplot















library("biomaRt")

# Select appropriate database
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)

# set up any filters and their appropriate values
filters <- c("hgnc_symbol")
gene <- c(A.1_list)
values <- as.list(gene)

# select attributes to return
attributes <- c("hgnc_id", "namespace_1003", "name_1006", "go_id", "go_linkage_type", "definition_1006")

# perform biomart query
bio <- getBM(attributes=attributes, filters=filters, values=values, mart=ensembl)
View(bio)



#Or Use Human Protein Atlas


install.packages("BiocManager")
## install hpar
BiocManager::install("hpar")
library("hpar")
allHparData()
data(hpaNormalTissue)

#Extract Column
data1 <- column$A.1
View(data1)

#Extract Column different way 
col <- df$A.1
col

col2 <- gconvert(col)
View(col2)

col3 <- col2$target
col3

col4 <- col2$target
newHPA <- getHpa(col4, hpadata= "rnaGeneTissue")
newHPA

#plot
plot <- barplot(newHPA$TPM, names.arg = newHPA$Tissue)

gg <- ggplot(newHPA, aes(x= Tissue, y= TPM)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))

gg

#Convert to Entrez
convert <- gconvert(data1)
View(convert)
class(convert)
target<- convert$target
class(target)
View(target)

#Get HPA data
HPA <- getHpa(target, hpadata= "hpaNormalTissue16.1")
HPA

class(HPA)


tissue <- table(HPA$Tissue)
tissue

class(tissue)

tissuedf <- as.data.frame(tissue)
tissuedf

hist <- barplot(tissuedf, height = tissuedf$Freq)

column$A.1

ggplot(tissuedf, aes(y=factor(tissuedf$Freq))) +geom_bar()
#image <- getHpa(target, type= "details")


#HPAStainR

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("HPAStainR")

library(HPAStainR)



HPA_data <- HPA_data_downloader(tissue_type = "both", save_file = TRUE)
HPA_data


listA.1 <- df$A.1
listA.1

listA.2 <- df$A.2
listA.2

stain <- HPAStainR::HPAStainR(listA.1, hpa_dat = HPA_data$hpa_dat)
stain

stain_df <- as.data.frame(stain)
stain_df
View(stain_df)

newstain_df <- stain_df[-c(122:184),]
View(newstain_df)
#Remove Rows

A.1plot <- ggplot(newstain_df, aes(x= cell_type, y= percent_high_expression)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1)) + xlab("A.1 Cluster")

A.1plot

finalstain_df <- stain_df[-c(41:184),]
finalstain_df

A.1plot1 <- ggplot(finalstain_df, aes(x= cell_type, y= percent_high_expression)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("A.1 Cluster High Expression")

A.1plot1


#A.2 
stain2 <- HPAStainR::HPAStainR(listA.2, hpa_dat = HPA_data$hpa_dat)
View(stain2)

stain2_df <- as.data.frame(stain2)
head(stain2_df)

#Search for all liver components and find very low expression

library(data.table)

table <- stain2_df[stain2_df$cell_type %like% "LIVER",]
View(table)

#Subset Rows

newstain2 <- stain2_df[-c(40:184),]
View(newstain2)

#GGplot to quantify A.2
View(newstain2)
A.2plot2 <- ggplot(newstain2,aes(x= cell_type, y= percent_high_expression)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("A.2 Cluster High Expression")

A.2plot2
