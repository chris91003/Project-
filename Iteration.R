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

class(column)

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



#Dataframe Exploration












#Tabular format of GO IDs and Intersecting Genes
intersect<- data.frame("Go_ID"= c(results$term_id), 
                        "Intersecting Genes"= c(results$intersection),
                       "Description" = c(results$term_name))
                       
intersect

View(intersect)






#Generate Manhattan Plots 

Ftable <- publish_gosttable(GO_F, use_colors = TRUE)
Ftable

GO_F$intersection

