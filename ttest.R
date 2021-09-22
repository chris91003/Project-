#T test to compare disease groups in phenotype list
library(ggplot2)
library(dplyr)

#read data
zscore <- read.csv("kendall-liver-cluster-profiles.csv")
zscore

phenotype <- read.csv("gtex_phenotypes_v8.csv")
View(phenotype)
View(zscore)
mean(zscore$A.1)

#Convert sample IDS into SUBJIDS
GTEx_SAMPID_to_SUBJID <- function(sampids){
  stringr::str_extract(string = sampids, pattern = "GTEX-[[:alnum:]]*")
}


zscore_SUBJID <- GTEx_SAMPID_to_SUBJID(zscore$X)
zscore_SUBJID



#Replace subject ID with new Subj ID and join
zscore$X <- zscore_SUBJID
View(zscore)

#rename zscore data column
names(zscore)[1] <- "SUBJID"

zscore

#Join dataframes together 
join <- left_join(zscore, phenotype, by= "SUBJID")

View(join)



#Explore Data




#Sort by Column A.1 by greater to less than

new_join <- join[order(-join$A.1),]
View(new_join)

order <- order(join[, "A.1"])
order
View(order)

#Extract just A.1 and hypertension columns and filter out 99 
mod_df <- data.frame("SUBJID"= new_join$SUBJID, 
                     "ZScore"= new_join$A.1, "Hypertension"= new_join$MHHTN)
mod_df
class(mod_df)

final_df <- subset(mod_df, mod_df$Hypertension != "99")
final_df

View(final_df)


#Visualize data

boxplot(final_df$ZScore~final_df$Hypertension)

barplot <- ggplot(final_df, aes(x= Hypertension))+
  geom_bar(fill= "orange") + xlab("Barplot of People with and without Hypertension")+
  ylab("Number of People")



#Subset top 25% and bottom 25% of data

top <- final_df[1:57, 1:3]
View(top)

hist(top$ZScore, main= "Zscore Distribution", xlab= "ZScores")


bottom <- final_df[169:226 ,1:3]
hist(bottom$ZScore, main= "Bottom Zscore Distribution", xlab= "ZScores")
View(bottom)



#Distribution is skewed so I used parametric tests
#Generate T test 
#Ho = mean of top zscore/ hypertension population is equal to mean of zscore and hypertension
#for bottom population

test_top <- t.test(top$ZScore, top$Hypertension)
test_top


test_bottom <- t.test(bottom$ZScore, bottom$Hypertension)
test_bottom

# 95 percent confident that means lie between 0.256 and 0.5313890


#Wilcox Test (assume that distribution isn't normal)

wilcox_top <- wilcox.test(top$ZScore, top$Hypertension)
wilcox_top

wilcox_bottom <- wilcox.test(bottom$ZScore, bottom$Hypertension)

hypertension_list <- c(test_top$p.value, test_bottom$p.value)

wilcox_table <- table(wilcox_top$p.value, wilcox_bottom$p.value)
wilcox_table

View(table)

column <- c(test_top$p.value, test_bottom$p.value)

column




#Liver Disease Analysis
join$MHLVRDIS




liver_disease <- data.frame("SUBJID"=join$SUBJID, 
                     "ZScore"= join$A.1, "Liver Disease"= join$MHLVRDIS)
View(liver_disease)

liver_disease_plot <- ggplot(liver_disease, aes(x= Liver.Disease))+
  geom_bar(fill= "orange") + xlab("Barplot of People with and without Liver Disease")+
  ylab("Number of People")

liver_disease_plot

top_liver <- liver_disease[1:57, 1:3]

bottom_liver <-  liver_disease[169:226 ,1:3]



liverdisease_test_top <- t.test(top_liver$ZScore, top_liver$Liver.Disease)
liverdisease_test_top

liverdisease_test_bottom <- t.test(bottom_liver$ZScore, bottom_liver
                                   $Hypertension)
liverdisease_test_bottom

liverdisease_list <- c(liverdisease_test_top$p.value, liverdisease_test_bottom$p.value)
liverdisease_list


Pvalue_df <- data.frame(unlist(hypertension_list), unlist(liverdisease_list))

names(Pvalue_df) <- c("P-value of Zscore/Hypertension", "P-value of 
                      Zscore/Liver Disease")
row.names(Pvalue_df)[1] <- "Top"
row.names(Pvalue_df)[2] <- "Bottom"
View(Pvalue_df)
