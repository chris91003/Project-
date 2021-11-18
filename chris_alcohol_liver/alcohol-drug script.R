library(tidyverse)
library(DESeq2)
```
setwd("/Users/cdima/OneDrive/Desktop/Project/chris_alcohol_liver")
Loading required data

## Required for GenesFromTissues
key_dat <- read.csv("data_in/key_table.csv")
## Has the phenotype data
phen_dat <- read.csv("data_in/gtex_phenotypes_v8.csv")

For more info on phenotypes go to the following website:
  https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000424/phs000424.v8.p2/pheno_variable_summaries/phs000424.v8.pht002742.v8.GTEx_Subject_Phenotypes.data_dict.xml


Loading personal functions


NewZScoreMaker <- function(gene_from_tissue_long_out){
  ## This function creates Z-scores based on the genes created from the long
  ## version of GenesFromTissue
  numeric_ind <- sapply(gene_from_tissue_long_out, is.numeric)
  zscores = rep(NA, nrow(gene_from_tissue_long_out))
  for (cur_tiss in unique(gene_from_tissue_long_out$tissue)) {
    tiss_ind <- gene_from_tissue_long_out$tissue %in% cur_tiss
    numeric_matrix <- gene_from_tissue_long_out[tiss_ind,numeric_ind]
    centered_matrix <- sweep(numeric_matrix, colMeans(numeric_matrix), MARGIN = 2)
    zscore_matrix <- sweep(centered_matrix,
                           apply(numeric_matrix, 2, sd),
                           FUN = "/",
                           MARGIN = 2)
    zscores[tiss_ind] <- apply(zscore_matrix, 1, sum, na.rm =T)/ncol(zscore_matrix)
  }
  return(zscores)    
  
}


GenesFromTissues <- function(tiss_of_int, genes_of_int, key_dat, vst_location,
                             run_all_available = FALSE){
  ### The purpose of this function is to easily load in normalized RNA-seq
  ### data to compare gene expression. The main idea is to see if the
  ### expression correlates across samples.
  ## load required packages
  require(DESeq2)
  require(tidyverse)
  ## Filter tissue based on any of the columns, this allows us to use abbreviations
  ## r names or names proper
  if (run_all_available == TRUE){
    tiss_of_int <- c()
    for (tissue in unique(key_dat$r_names)) {
      tiss_path <- paste0(vst_location,tissue,"-vsd-mean-filtered.rda")
      if (file.exists(tiss_path)) {
        tiss_of_int <- append(tiss_of_int, tissue)
      }
    }
    
  }
  cur_keys <- key_dat %>% select(abrreviated, r_names, official_name) %>%
    unique() %>%
    filter(case_when(
      abrreviated %in% tiss_of_int |
        r_names %in% tiss_of_int |
        official_name %in% tiss_of_int ~ TRUE,
      TRUE ~ FALSE
    ))
  cur_keys <- cur_keys[order(tiss_of_int, cur_keys$r_names),]
  
  ## Prepare returned data frame
  final_gene_dat <- NULL
  for (tiss_ind in seq_along(tiss_of_int)) {
    ## Select and load tissues which exist
    tissue <- cur_keys$r_names[tiss_ind]
    abbrev <- cur_keys$abrreviated[tiss_ind]
    tiss_path <- paste0(vst_location,tissue,"-vsd-mean-filtered.rda")
    ## Catch missing vsts
    if (!file.exists(tiss_path)) {
      stop(paste0(tissue, " PATH does not exist"))
    }
    ## Load data and do a quick sanity check
    load(tiss_path)
    if(!all(gtabMeanFiltered$gene_id == rownames(assay(generalVSDMeanFiltered)))){
      stop(paste0(tissue, " fails sanity check"))
    }
    ## Index on genes of interest
    ind <- gtabMeanFiltered$gene_name %in% genes_of_int
    final_gtab <- gtabMeanFiltered[ind,]
    ## Take out gene data and name it
    gene_dat <- t(assay(generalVSDMeanFiltered)[ind,,drop =FALSE])
    colnames(gene_dat) <- final_gtab$gene_name
    ## If more than one tissue bind the tissue name
    if (length(tiss_of_int) > 1) {
      gene_dat <- cbind(as.data.frame(gene_dat), tissue)
      ## if first tissue just place it in final object
      if (is.null(final_gene_dat)) {
        final_gene_dat <- gene_dat
        next
      }
      #Use bind rows from dplyr in case not all genes are shareed
      final_gene_dat <- bind_rows(final_gene_dat, gene_dat) %>%
        select(everything(), tissue)
    } else{
      colnames(gene_dat) <-final_gtab$gene_name
      gene_dat <- rownames_to_column(as.data.frame(gene_dat), var = "SAMPID")
      final_gene_dat <- gene_dat
    }
    
  }
  if (!("SAMPID" %in% colnames(final_gene_dat))) {
    final_gene_dat <- final_gene_dat %>% rownames_to_column(var = "SAMPID")
  }
  ## 
  final_gene_dat <- final_gene_dat  %>%
    mutate(SUBJID = GTEx_SAMPID_to_SUBJID(SAMPID)) %>%
    select(SUBJID, SAMPID, everything())
  
  
  
  return(final_gene_dat)
  
}
## Makes SAMPID names SUBJID ma,es
GTEx_SAMPID_to_SUBJID <- function(sampids){
  ## A simple function to quickly turn SAMPIDs into SUBJIDs to connect
  ## individuals to their phenotype
  stringr::str_extract(string = sampids, pattern = "GTEX-[[:alnum:]]*")
}

GeneClusterPuller <- function(r_name_vect, location){
  ### This function is used to pull the genes of tissue clusters
  tiss_list <- list()
  for (tiss in r_name_vect) {
    cur_frame <- read.csv(paste0(location,
                                 "kendall-",
                                 tiss,
                                 "-gene-clusters-high-variance.csv"))
    clust_list <- list()
    for (cur_col in seq(ncol(cur_frame))) {
      gene_column <- cur_frame[,cur_col]
      gene_column <- gene_column[gene_column !=""]
      column_name <- colnames(cur_frame)[cur_col]
      clust_list[[column_name]] <- gene_column
    }
    tiss_list[[tiss]] <- clust_list
  }
  return(tiss_list)
}

#### Load in liver A.1 cluster genes

liver_clusters <- GeneClusterPuller("liver","data_in/variance_genes_kendall/")
a1_cluster <- liver_clusters$liver$A.1

a1_cluster_df <- as.data.frame(a1_cluster)
a1_cluster_df

#Genes of Interest
library(tidyverse)
CYP_genes <- startsWith(a1_cluster,'CYP')
CYP_genes

CYP <- a1_cluster_df %>% filter(str_detect(a1_cluster, "CYP"))

ADH_genes <- startsWith(a1_cluster, 'ADH')

ADH <- a1_cluster_df %>% filter(str_detect(a1_cluster, "ADH"))


genes_of_interest_CYP <- c(CYP)

genes_of_interest_ADH <- c(ADH)

genes_of_interest <- c("ADH4", "ADH1A", "ADH1C")

## Load data
liver_data <- GenesFromTissues("liver", genes_of_interest, key_dat, "data_in/",
                               run_all_available = F)


liver_data$tissue <- "liver"

## This will make a Zscore from all of the genes you loaded in above
liver_data$zscore <- NewZScoreMaker(liver_data)


liver_data$zscore

liver_data

### Now join the data with the phenotype data and use your skills to find stuff out!




#Joining phenotype data with ADH genes

zscore <- liver_data$zscore
zscore

phenotype <- read.csv("gtex_phenotypes_v8.csv")
View(phenotype)
View(zscore)
mean(zscore$A.1)



#Join dataframes together (ADH data)


ADH_data <- left_join(liver_data, phenotype, by= "SUBJID")
ADH_data$DTHHRDY <- as.factor(ADH_data$DTHHRDY)

hardy_factor <- factor(ADH_data$DTHHRDY, levels = c("1", "2", "3", "4", "0"))

View(ADH_data)

#ADH data viszualization

adh_plot <- ggplot(ADH_data, aes(x=DTHHRDY, y=zscore))+
  geom_violin(aes(fill = DTHHRDY))+geom_boxplot(width=.3)+geom_jitter()+ labs(y
                                                                              = "ADH Zscore")


adh_jitter <- ggplot(ADH_data) + aes(DTHHRDY, zscore) +
  geom_jitter()+labs(y= "ADH Zscore")

#Model
ADH_model <- lm(ADH_data$zscore~ADH_data$DTHHRDY)

summary(ADH_model)


#Deathplace Visualization
ADH_data$DTHPLCE

death <- ggplot(ADH_data, aes(x= DTHPLCE, fill= DTHPLCE))+
  geom_bar()+ theme_gray()

death+coord_flip() +labs(x= "Deathplace", y= "Number of People")



##Filter by place of death DTHPLCE

ADH_data$DTHPLCE == "Hospital inpatient"
hospital <- ADH_data[grep("Hospital inpatient", ADH_data$DTHPLCE), ]
hospital$DTHPLCE

View(hospital)

hospital_plot <- ggplot(hospital) + aes(DTHHRDY, zscore) +
  geom_jitter()

summary(lm(hospital$zscore~hospital$DTHHRDY))


#Deathplace outside Hospital

home <- ADH_data[grep("Decedent's home", ADH_data$DTHPLCE), ]
View(home)

summary(lm(home$zscore~home$DTHHRDY))


home_plot <- ggplot(home) + aes(DTHHRDY, zscore) +
  geom_jitter()


home1 <- ADH_data[grep("Decedent's home", "Emergency room", ADH_data$DTHPLCE), ]


#Correlate Drinking with HARDY
ADH_data$MHDRNKSTS <- ifelse(ADH_data$MHDRNKSTS == "Yes", 1, 0)

ADH_data$MHDRNKSTS <- as.factor(ADH_data$MHDRNKSTS)

ggplot(ADH_data) + aes(MHDRNKSTS, zscore) +
  geom_jitter()


#Plotting zscore vs alcohol use for people in hospital

hospital$MHDRNKSTS <- ifelse(hospital$MHDRNKSTS == "Yes", 1, 0)

hospital$MHDRNKSTS <- as.factor(hospital$MHDRNKSTS)


ggplot(hospital) + aes(MHDRNKSTS, zscore) +
  geom_jitter()

View(hospital)



#Plotting zscore vs alcohol use for people not in hospital

ADH_data$DTHPLCE == "Emergency room"


non_hosp <-ADH_data[(ADH_data$DTHPLCE == "Emergency room") 
                    | (ADH_data$DTHPLCE == "Decendent's home"),] 

ggplot(non_hosp, aes( MHDRNKSTS, zscore)) +geom_jitter() +labs(x= "Nonhospital Drinkers")


#Joining phenotype data with CYP genes
CYPgenes <- c("CYP4A11", "CYP4A22", "CYP3A7", "CYP3A4", "CYP3A43",
              "CYP7A1", "CYP2C19", "CYP2C8", "CYP1A1", "CYP1A2",
              "CYP2A6", "CYP2A7", "CYP2B7P", "CYP2B6", "CYP2A13")

liver_data2 <- GenesFromTissues("liver", CYPgenes, key_dat, "data_in/",
                               run_all_available = F)


liver_data2$tissue <- "liver"

## This will make a Zscore from all of the genes you loaded in above
liver_data2$zscore <- NewZScoreMaker(liver_data2)


liver_data2$zscore



#Joining phenotype data with ADH genes

zscore2 <- liver_data2$zscore
zscore2

phenotype <- read.csv("gtex_phenotypes_v8.csv")
View(phenotype)
View(zscore)
mean(zscore$A.1)


#Join dataframes together 


CYP_data <- left_join(liver_data2, phenotype, by= "SUBJID")
CYP_data$DTHHRDY <- as.factor(CYP_data$DTHHRDY)

View(CYP_data)


#Plot CYP
cyp_plot <- ggplot(CYP_data, aes(x=DTHHRDY, y=zscore))+
  geom_violin(aes(fill = DTHHRDY))+geom_boxplot(width=.3)+geom_jitter() +labs(y= "CYP Zscore")

cyp_jitter <- ggplot(ADH_data) + aes(DTHHRDY, zscore) +
  geom_jitter()+labs(y= "CYP Zscore")

#CYP model

CYP_model <- lm(CYP_data$zscore~CYP_data$DTHHRDY)
summary(CYP_model)



library(gridExtra)
grid.arrange(adh_plot, cyp_plot, nrow=1)
grid.arrange(adh_jitter, cyp_jitter, nrow=1)                    
