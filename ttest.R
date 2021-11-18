#T test to compare disease groups in phenotype list
library(ggplot2)
library(dplyr)
library(ggthemes)
setwd("/Users/cdima/OneDrive/Desktop/Project")
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


#Trying to obtain all columns with sum of >1 = 100 or more
empty_list <- c()

for (i in ncol(join)){
  empty_list[i] <- sum(join[, i] == 1)
  
}

empty_list

library(dplyr)
join %>%
summarise_if(is.numeric, sum, na.rm= TRUE)





#Explore Data
join$



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

boxplot(join$ZScore~final_df$Hypertension)

barplot <- ggplot(join, aes(x= Hypertension))+
  geom_bar(fill= "blue") + xlab("Barplot of People with and without Hypertension")+
  ylab("Number of People")

barplot





#Distribution is skewed so I used parametric tests
#Generate T test 
#Ho = mean of top zscore/ hypertension population is equal to mean of zscore and hypertension
#for bottom population


#Ttest smokers
smokers <- ifelse(join$MHSMKSTS == "Yes", 1, 0)
smokers
smoke <- t.test(join$A.1, smokers)
smoke
smoke$p.value


t.test(join$A.1, smokers)$p.value

smoker_test <- t.test(join$A.1, smokers)
smoker_test
smoker_pval <- smoker_test$p.value
smoker_pval

#Hypertension
hyp_df <- subset(join, join$MHHTN != "99")
hyp_df

#Hypertension
hyp_df$MHHTN

hyp_test <- t.test(join$A.1, hyp_df$MHHTN)
hyp_test

hyp_pval <- hyp_test$p.value
hyp_pval

#Ttest Pancreas Contamination and Hypertension

contamination <- t.test(join$B,hyp_df$MHHTN)
contamination

contamination_pval <- contamination$p.value
contamination_pval

#table
ne <- table(hyp_pval, contamination_pval, smoker_pval)
new <- as.data.frame(ne)

View(new)
final <- subset(new, select= -c(Freq))
final
names(final) <- c("Hypertension/Drug Metabolizing Cluster", " Hypertension/Contamination Cluster", 
                  "Smoker/Drug Metabolizing Cluster")
final

row.names(final)[1] <- "P-Value"
final
View(final)

View(join)
hyp_df$MHHTN

#
A.2_pval <- t.test(join$A.2, hyp_df$MHHTN)$p.value
C_pval <- t.test(join$C, hyp_df$MHHTN)$p.value
D_pval <- t.test(join$D, hyp_df$MHHTN)$p.value
E_pval <- t.test(join$E, hyp_df$MHHTN)$p.value

hyp_df$MHHTN

y <- table(A.2_pval, C_pval, D_pval, E_pval)

as.data.frame(y)
View(as.data.frame(y))

hist(join$A.1)
hist(join$A.2)
hist(join$B)



#Wilcox Test (assume that distribution isn't normal)

wilcox_top <- wilcox.test(top$ZScore, top$Hypertension)
wilcox_top

wilcox_bottom <- wilcox.test(bottom$ZScore, bottom$Hypertension)

hypertension_list <- c(test_top$p.value, test_bottom$p.value)
hypertension_list


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




#T test showed correlations that shouldn't have worked. Maybe Regression Analysis
#is better model?

join$MHSMKSTS <- ifelse(join$MHSMKSTS == "Yes", 1, 0) 





#Basic Linear Model Between A.1 and Smokers

linear_smokers <- lm(join$MHSMKSTS~join$A.1, data= join)
linear_smokers

summary(linear_smokers)

rsquared <- summary(linear_smokers)$r.squared

R2 <- round(rsquared, digits =4)
R2

#Ggplot A.1 and Smoking  

smoke_A.1plot <- ggplot(join ,aes(x= A.1, y= MHSMKSTS), xlab= "Liver Drug Metabolizing Cluster",
       ylab= "Smokers vs NonSmokers") + geom_point() + geom_smooth(formula = y~x, method= "lm") +
  theme_base() + annotate("text", x= 0.8, y= 1.2, label= paste0("R-Squared: ", R2))

smoke_A.1plot+ labs(x= "A.1 Zscore", y= "Smoker vs NonSmoker")



summary(ln)


#Regular plot

plot(x= join$A.1, y= join$MHSMKSTS, xlab= "Cluster A.1 Zscores", ylab=
       "Smokers vs NonSmokers")

abline(linear, col= "red")



#Linear Model Hypertension

View(join)

linear_hypertension <- lm(join$MHHTN~join$A.1, data= join)


summary(linear_hypertension)

join$MHHTN <- replace(join$MHHTN, join$MHHTN == 99, NA)

join$MHHTN

join$MHHTN <- factor(join$MHHTN)
join$MHHTN

model <- lm(join$A.1~join$MHHTN)
summary(model)


boxplot(A.1~ MHHTN, data= join)

join$MHHTN

summary(new_model)

other <- lm(factor~join$A.1)


#Logistic Regression Analysis
summary(join)

#Convert to Factors

join$MHSMKSTS <- ifelse(join$MHSMKSTS == "Yes", 1, 0)

join$MHSMKSTS
join$A.1 <- as.factor(join$A.1)

#Visualize Data
table(join$A.1, join$MHSMKSTS)
plot(join$A.1, join$MHSMKSTS)
cor(join$A.1, join$MHSMKSTS)


#Simple Logistic Regression model

logreg <- glm(join$MHSMKSTS~join$A.1, family= "binomial", data= join)

logreg


summary(logreg)

#Partition into Training and Test Set
library(caret)

set.seed(1234)

train <- createDataPartition(join, p= .6, list= FALSE)



#Multi Linear Regression Smokers, Age, A.1

plot(join$A.1, join$AGE)


Smoker_Age_A.1 <- lm(join$A.1~ join$MHSMKSTS+ join$AGE, data= join)


summary(Smoker_Age_A.1)

#Multi Linear Regression Model Hypertension, Smokers, A.1, Age

join[join == 99] <- NA

join$MHHTN <- replace(join$MHHTN, join$MHHTN == 99, NA)

join$MHHTN

join1

join1[join1 == 99] <- NA


which(is.na(join$MHHTN), arr.ind= TRUE)
join1 <- join[!is.na(join$MHHTN), ]

join1_mod <- join1[-c(123,213),]

join1_mod$H


A.1_Smoker_MHHTN<- lm(join1_mod$A.1~join1_mod$MHHTN+join1_mod$MHSMKSTS, data= join1_mod)
summary(A.1_Smoker_MHHTN)




#Remove NA from dataframe

join[join == 99] <- NA

#A2 Cluster Regression with Post-Mortom Interval time? (time person died to time of autopsy) (TRISCHD)

hist(join$TRISCHD)



summary(lm(join$A.2~join$TRISCHD))


#Visualization of Ischemic Time, A.2 Cluster, and Hardy

ggplot(join, aes(TRISCHD, A.2))+
  geom_point()+geom_smooth(method= 'lm', se= FALSE)

ggplot(join, aes(TRISCHD, A.2, color= DTHHRDY))+
  geom_point()

ggplot(join, aes(TRISCHD, A.2))+
  geom_point()+
  facet_wrap(~DTHHRDY, nrow = 1)


cor(join$A.2, join$TRISCHD)


#A.2 Cluster, Ischemic time, Hardy Model
join$DTHHRDY <- as.factor(join$DTHHRDY)

A.2_model <- lm(A.2~TRISCHD +DTHHRDY, data= join)


#C Cluster and Ischemic Time

ggplot(join, aes(TRISCHD, C))+geom_jitter()+geom_smooth(method= 'lm', se= FALSE)

ggplot(join, aes(DTHHRDY, C)) +geom_jitter(width= 0.1)

C_model <- lm(C~TRISCHD +DTHHRDY, data= join)

cor(join$C,join$TRISCHD)



#E cluster Contamination


#F Cluster

join$SEX <- as.factor(join$SEX)
ggplot(join, aes(SEX, F))+geom_jitter(width= 0.1)


#G cluster and hypertension (blood pressure)

join$MHHTN <- as.factor(join$MHHTN)
ggplot(join, aes(MHHTN, G)) +
  geom_jitter(width= 0.1)


join$MHCVD <- as.factor(join$MHCVD)
ggplot(join, aes(	MHHTN, G)) + geom_jitter(width= 0.1)

#H Cluster maybe liver injury?




#I.2 cluster and sex
join$SEX <- as.factor(join$SEX)


ggplot(join, aes(SEX, I.2))+geom_jitter(width= 0.2, aes(color = SEX))





#Liver Cluster T2D

plot(join$A.1, join$MHT2D)
plot(join$MHT2D, join$A.1)


diabetes <- lm(join$MHT2D~ join$A.1)

summary(diabetes)
rsquared <- summary(diabetes)$r.squared


#Multi Regression Analysis

diabetes_age <- lm(join$MHT2D~join$A.1+join$AGE)
summary(diabetes_age)


diabetes_age_weight <- lm(MHT2D~A.1+AGE+WGHT, data= join)
summary(diabetes_age_weight)




summary(lm(join$A.1~join$MHT2D+join$AGE+join$WGHT))


#library(car)


avPlots(diabetes_age_weight)

#plot(diabetes_age_weight)




#HardyScale Analysis 0= ventilator case 1= fast/violent death 
#2- fast natural death 3- intermediate death 4-slow death

#Remove NA

hardy <- join[!is.na(join$DTHHRDY),]

boxplot(A.1~DTHHRDY, data= hardy)

#Violin plot Hardy Score

join$DTHHRDY <- as.factor(join$DTHHRDY)
join$DTHHRDY

hardy_plot <- ggplot(join, aes(x=DTHHRDY, y=A.1))+
  geom_violin(aes(fill = DTHHRDY))+geom_boxplot(width=.3)
  

hardy_plot+geom_jitter()


#Violin plot without the NA

View(hardy)
hardy$DTHHRDY <- as.factor(hardy$DTHHRDY)
hardy1_plot <- ggplot(hardy, aes(x=DTHHRDY, y=A.1))+
  geom_violin(aes(fill = DTHHRDY))+geom_boxplot(width=.2)+
  labs(title = "Hardy Score and A.1 Cluster Comparison",
       x= "Hardy Score", y= "Cluster Z-Score") +geom_jitter()


hardy1_plot+geom_jitter()


hardyA.2_plot <- ggplot(hardy, aes(x=DTHHRDY, y=A.2))+
  geom_violin(aes(fill = DTHHRDY))+geom_boxplot(width=.2)+
  labs(title = "Hardy Score and A.2 Cluster Comparison",
       x= "Hardy Score", y= "Cluster Z-Score")+geom_jitter()


hardyA.2_plot+geom_jitter()



library(gridExtra)
grid.arrange(hardy1_plot, hardyA.2_plot, nrow= 1)
#Hardy Model A1 and A2

join$DTHHRDY <- factor(join$DTHHRDY, levels = (1, 2, 3, 4, 0))
join$DTHHRDY
summary(lm(A.1~DTHHRDY, data= join))

hardy_model <- lm(A.1~DTHHRDY, data= hardy)
summary(hardy_model)


plot(density(resid(hardy_model))
qqnorm(resid(hardy_model))


hardy_modelA.2 <- lm(A.2~DTHHRDY,data= hardy)
summary(hardy_modelA.2)     

hardy_model$coefficients


cor(join$DTHHRDY, join$A.1, use = "complete.obs")



#Hardy Model C

C_hardy <- ggplot(join) + aes(DTHHRDY, C) +
  geom_jitter()

D_hardy <- ggplot(join) + aes(DTHHRDY, D) +
  geom_jitter()

#Diabetes
ggplot(join) + aes(MHT2D, A.2) +
  geom_jitter()

