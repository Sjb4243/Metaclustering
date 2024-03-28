library(GEOquery)
library(caret)
library(glmnet)
library(tidyverse)
library(factoextra)

setwd('/Users/u5590987/Documents/Animesh/metaclustering')

#This script does require some user input namely: 
#Many of the variable names are for the GSE25097 and GSE63898
#A desired number of clusters is necessary for the clustering stage
#Once the data has been prepped, the rest of the pipeline is automatic

#DATA PREP
#Need to delineate between healthy/tumor samples
######################################################################################

GSE63898 <- getGEO("GSE63898", GSEMatrix =TRUE, getGPL=FALSE)
if (length(GSE63898) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
GSE63898 <- GSE63898[[idx]]
GSE63898_samples <- GSE63898$title

GSE25097 <- getGEO("GSE25097", GSEMatrix =TRUE, getGPL=FALSE)
if (length(GSE25097) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
GSE25097 <- GSE25097[[idx]]
GSE25097_samples <- GSE25097$title

GSE25097 <- read.delim('data/GSE25097.gene.txt')
GSE63898 <- read.delim('data/GSE63898.gene.txt')
######################################################################################

#Filtering to the tumour samples for both datasets
GSE25097_tumour <- GSE25097[,-c(1:2)]
GSE25097_tumour <- GSE25097_tumour[,startsWith(GSE25097_samples, 'tumor')]
GSE25097_tumour <- cbind(Symbol = GSE25097[,2], GSE25097_tumour)

GSE63898_tumour <- GSE63898[,-c(1:2)]
GSE63898_tumour <- GSE63898_tumour[,startsWith(GSE63898_samples, 'hepato')]
GSE63898_tumour <- cbind(Symbol = GSE63898[,2], GSE63898_tumour)


######################################################################################
#Processing the data as described in the paper
#HAve to pivot to long, then cast to wide to get samples to be the rows and genes the columns
#I would convert this into a function for all datasets but the instructions are slightly different for each
#Near zero var didnt work for GSE63898, discussed in paper
processed_GSE25097 <- GSE25097_tumour %>% filter(Symbol != '---')
processed_GSE25097<- processed_GSE25097%>% pivot_longer(cols = 2:ncol(GSE25097_tumour), names_to = 'sample', values_to = 'value')
processed_GSE25097 <- processed_GSE25097 %>% pivot_wider(names_from = Symbol, values_from = value)
GSE25097_patients <- processed_GSE25097$sample
processed_GSE25097 <- round(processed_GSE25097[,-1],1)
processed_GSE25097 <- processed_GSE25097[,-(nearZeroVar(processed_GSE25097))]
processed_GSE25097 <- scale(processed_GSE25097[,-1])
processed_GSE25097 <- as.data.frame(processed_GSE25097)
rownames(processed_GSE25097) <- GSE25097_patients


processed_GSE63898 <- GSE63898_tumour %>% filter(Symbol != '---')
processed_GSE63898 <- processed_GSE63898 %>% filter(Symbol != 'SLC35E2')
processed_GSE63898<- processed_GSE63898%>% pivot_longer(cols = 2:ncol(GSE63898_tumour), names_to = 'sample', values_to = 'value')
processed_GSE63898 <- processed_GSE63898 %>% pivot_wider(names_from = Symbol, values_from = value)
GSE63898_patients <- processed_GSE63898$sample
processed_GSE63898 <- round(processed_GSE63898[,-1],1)
processed_GSE63898 <- scale(processed_GSE63898[,-1])
processed_GSE63898 <- as.data.frame(processed_GSE63898)
rownames(processed_GSE63898) <- GSE63898_patients

############################################################################################
#At this stage, the data consists solely of genes in the columns and is ready for clustering
############################################################################################
#Set cluster number based on WCSS plot 
cluster_number = 3
GSE25097.model = kmeans(processed_GSE25097, cluster_number, nstart = 20)
GSE63898.model = kmeans(processed_GSE63898, cluster_number, nstart = 20)

############################################################################################
                                  #Cluster visualisation
############################################################################################

plot1 <- fviz_cluster(GSE25097.model, data = processed_GSE25097, geom = 'point') + 
  labs(title = 'GSE25097 clustering of processed data') + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave('GSE25097clusters.png')

plot2 <- fviz_cluster(GSE63898.model, data = processed_GSE63898, geom = 'point') + 
  labs(title = 'GSE63898 clustering of processed data') + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave('GSE25097clusters.png')

############################################################################################
                                  #Elastic net
############################################################################################
#CAUTION - I've found that running this in rstudio is impossible due to memory
#Run through command line instead

#Adding cluster data to original datasets
GSE25097_elastic <- cbind(processed_GSE25097, cluster_id = GSE25097.model$cluster)
GSE63898_elastic <- cbind(processed_GSE63898, cluster_id = GSE63898.model$cluster)

#Function to perform elastic net with whatever dataset you give it
#Make sure the cluster_id column is called cluster_id
perform_elastic <- function(dataset){
  dataset$cluster_id <- as.factor(dataset$cluster_id)
  indexes = createDataPartition(dataset$cluster_id, p = 0.75, list = FALSE)
  trn = dataset[indexes, ]
  tst = dataset[-indexes, ]
  cv_5 = trainControl(method = "cv", number = 5)
  elnet = train(
    cluster_id ~ ., data = trn,
    method = "glmnet",
    trControl = cv_5
  )
  return_list <- list(elnet, tst)
  return(return_list)
}
#Function returns both the model(1st) and the test data for testing (2nd)
#Test data isnt used here, but can be used to determine accuracy

returned_list_GSE25097 <- perform_elastic(GSE25097_elastic)
model_GSE25097 <- returned_list_GSE25097[[1]]
test_data_GSE25097 <- returned_list_GSE25097[[2]]
importances <- varImp(model_GSE25097, lamba = model_GSE25097$lamba.min)
GSE25097_importances <- importances$importance
GSE25097_importances$genes <- rownames(GSE25097_importances)

returned_list_GSE63898 <- perform_elastic(GSE63898_elastic)
model_GSE63898 <- returned_list_GSE63898[[1]]
test_data_GSE63898 <- returned_list_GSE63898[[2]]
importances <- varImp(model_GSE63898, lamba = model_GSE63898$lamba.min)
GSE63898_importances <- importances$importance
GSE63898_importances$genes <- rownames(GSE63898_importances)


getElasticResults <- function(dataset, topx){
  datasetName <- as.character(substitute(dataset))
  for (cluster in 1:cluster_number){
    sorted_data <- dataset[order(-dataset[,cluster]),]
    top_genes <- paste(sorted_data[1:topx,4], collapse = "\n")
    cat('\n', 'Top genes for cluster', cluster, '\n')
    cat(top_genes)
    write(top_genes, file = paste(datasetName, 'cluster', cluster, '.txt'))
  }
}
getElasticResults(GSE25097_importances,30)
getElasticResults(GSE63898_importances,30)


