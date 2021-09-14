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



HPA_data <- HPA_data_downloader(tissue_type = "both", save_file = FALSE)
HPA_data


listA.1 <- df$A.1
listA.1

listA.2 <- df$A.2
listA.2

stain <- HPAStainR::HPAStainR(listA.1, hpa_dat = HPA_data$hpa_dat)
View(stain)

stain_df <- as.data.frame(stain)
stain_df
View(stain_df)

newstain_df <- stain_df[-c(122:184),]
View(newstain_df)











#A.1 Plot

A.1plot <- ggplot(newstain_df, aes(x= cell_type, y= percent_high_expression)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1)) + xlab("A.1 Cluster")

A.1plot

#Remove Rows
finalstain_df <- stain_df[-c(41:184),]
finalstain_df

A.1plot1 <- ggplot(finalstain_df, aes(x= cell_type, y= percent_high_expression, fill= cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("A.1 Cluster High Expression") 
  

A.1plot1


#A.2 
stain2 <- HPAStainR::HPAStainR(listA.2, hpa_dat = HPA_data$hpa_dat)
View(stain2)

stain2_df <- as.data.frame(stain2)
head(stain2_df)

#Search for all liver components and find very low expression

library(data.table)

A.2table <- stain2_df[stain2_df$cell_type %like% "LIVER",]
View(A.2table)

#Subset Rows

newstain2 <- stain2_df[-c(40:184),]
View(newstain2)

#GGplot to quantify A.2
View(newstain2)
A.2plot2 <- ggplot(newstain2,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("A.2 Cluster High Expression")


  A.2plot2  


  


  
  
  
  
#B stain
Bstain <- HPAStainR::HPAStainR(df$B, hpa_dat = HPA_data$hpa_dat)
View(Bstain)  
  
B_stain <- as.data.frame(Bstain)

Bmod <- B_stain[-c(11:184),]
Bmod  
  
  
Bplot <- ggplot(Bmod,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("B Cluster High Expression")  
  
  
Bplot  

#Bplot Liver expression
Btable <- B_stain[B_stain$cell_type %like% "LIVER",]
View(Btable)  
  




#C stain

Cstain <- HPAStainR::HPAStainR(df$C, hpa_dat = HPA_data$hpa_dat)
View(Cstain)  

C_stain <- as.data.frame(Cstain)

Cmod <- C_stain[-c(16:184),]
Cmod  


Cplot <- ggplot(Cmod,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("C Cluster High Expression")  


Cplot  

#Cplot Liver expression
Ctable <- C_stain[C_stain$cell_type %like% "LIVER",]
View(Ctable)  




#D stain

Dstain <- HPAStainR::HPAStainR(df$D, hpa_dat = HPA_data$hpa_dat)
View(Dstain)  

D_stain <- Dstain[-c(75:184),]

Dplot <- ggplot(D_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("D Cluster High Expression")  


Dplot  

#Dplot Liver expression
Dtable <-D_stain[D_stain$cell_type %like% "LIVER",]
View(Dtable)  



#E stain

Estain <- HPAStainR::HPAStainR(df$E, hpa_dat = HPA_data$hpa_dat)
View(Estain)  

E_stain <- Estain[-c(11:184),]

Eplot <- ggplot(E_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("E Cluster High Expression")  

Eplot


#F stain


Fstain <- HPAStainR::HPAStainR(df$F, hpa_dat = HPA_data$hpa_dat)
View(Fstain)  

F_stain <- Fstain[-c(41:184),]

Fplot <- ggplot(F_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("F Cluster High Expression")  



Fplot



#G stain

Gstain <- HPAStainR::HPAStainR(df$G, hpa_dat = HPA_data$hpa_dat)
View(Gstain)  

G_stain <- Gstain[-c(26:184),]

Gplot <- ggplot(G_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("G Cluster High Expression")  

Gplot


#G plot for low expression 

Glow <- ggplot(G_stain,aes(x= cell_type, y= percent_low_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("G Cluster Low Expression")  


Glow



#H stain

Hstain <- HPAStainR::HPAStainR(df$H, hpa_dat = HPA_data$hpa_dat)
View(Hstain)  

H_stain <- Hstain[-c(21:184),]

Hplot <- ggplot(H_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("H Cluster High Expression")  

Hplot


#I.1 stain

I.1stain <- HPAStainR::HPAStainR(df$I.1, hpa_dat = HPA_data$hpa_dat)
View(I.1stain)  

I.1_stain <- I.1stain[-c(21:184),]

I.1plot <- ggplot(I.1_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("I.1 Cluster High Expression")  

I.1plot


#I.2 stain


I.2stain <- HPAStainR::HPAStainR(df$I.2, hpa_dat = HPA_data$hpa_dat)
View(I.2stain)  

I.2_stain <- I.2stain[-c(31:184),]

I.2plot <- ggplot(I.2_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("I.2 Cluster High Expression")  

I.2plot


#I.2 plot low

I.2low_stain <- I.2stain[-c(55:184),]

I.2low <- ggplot(I.2low_stain,aes(x= cell_type, y= percent_low_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("I.2 Cluster Low Expression")  


I.2low











#Attach all graphs


par(mfrow=c(5,6))
a1 <- ggplot(finalstain_df, aes(x= cell_type, y= percent_high_expression, fill= cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("A.1 Cluster High Expression") 

newa1 <- a1 +theme(legend.position = "none")

newa1
a2 <- ggplot(newstain2,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("A.2 Cluster High Expression")

newa2 <- a2 +theme(legend.position = "none")

b<- ggplot(Bmod,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("B Cluster High Expression")  

newb <- b +theme(legend.position = "none")

c <- ggplot(Cmod,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("C Cluster High Expression")  

newc <- c +theme(legend.position = "none")

d<- ggplot(D_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("D Cluster High Expression")  

newd <- d +theme(legend.position = "none")

e<- ggplot(E_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("E Cluster High Expression")  

newe <- e +theme(legend.position = "none")

f<- ggplot(F_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("F Cluster High Expression")  

newf <- f +theme(legend.position = "none")

g<- ggplot(G_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("G Cluster High Expression")  

newg <- g +theme(legend.position = "none")

h<-  ggplot(H_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("H Cluster High Expression")  

newh <- h +theme(legend.position = "none")

i.1 <- ggplot(I.1_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("I.1 Cluster High Expression")  

newi.1 <- i.1 +theme(legend.position = "none")


i.2 <-ggplot(I.2_stain,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("I.2 Cluster High Expression")  

newi.2 <- i.2 +theme(legend.position = "none")




d1 <- Dstain[-c(61:184),]
d2<- ggplot(d1,aes(x= cell_type, y= percent_high_expression, fill=cell_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))+ xlab("D Cluster High Expression")  

d2
newd2 <- d2 +theme(legend.position = "none")


install.packages("gridExtra")
library(gridExtra)
library(ggplot2)
grid <- grid.arrange(newa1, newa2, newb, newc, newd, newe, newf , newg,
             newh, newi.1, newi.2, nrow= 5, ncol= 6, widths= c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1))

grid1 <- grid.arrange(newa1, newa2, newb, newc, newd, newe, widths= c(1, 1, 1))

grid2 <- grid.arrange(newf, newg, newh, newi.2, nrow= 1, ncol= 4, widths= c(1,1,1,1))
View(grid)


newgrid <- grid.arrange(newa1, newa2, newb, newc, newd2, newe, widths= c(0.6, 0.6, 0.6))





#Loop to perform stain for every cluster    
for i in 1:length(split1)) {
  genstain <- HPAStainR::HPAStainR(split1[[i]], hpa_dat = HPA_data$hpa_dat)
  print(genstain)
}  

View(genstain)


#Maybe lapply
lapply(df, function(x, HPAStainR::HPAStainR(split1[[i]], hpa_dat = HPA_data$hpa_dat) )

  