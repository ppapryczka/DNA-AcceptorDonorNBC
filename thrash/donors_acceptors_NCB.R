#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# ---------------------------------------------------------------------------------------------
# LOAD PACKAGES 
# ---------------------------------------------------------------------------------------------

# stringr - for some string functions
if (! "stringr" %in% row.names(installed.packages()))
  install.packages("stringr")
library(stringr)

# purrr - map function
if (! "purrr" %in% row.names(installed.packages()))
  install.packages("purrr")
library(purrr)

# e1071 - naive bayes
if (! "e1071" %in% row.names(installed.packages()))
  install.packages("e1071")
library(e1071)

# caret - confusion matrix
if (! "caret" %in% row.names(installed.packages()))
  install.packages("caret")
library(caret)

if (! "randomForest" %in% row.names(installed.packages()))
  install.packages("randomForest")
library(randomForest)

if (! "rsample" %in% row.names(installed.packages()))
  install.packages("rsample")
library(rsample)

if (! "Hmisc" %in% row.names(installed.packages()))
  install.packages("Hmisc")
library(Hmisc)

if (! "corrplot" %in% row.names(installed.packages()))
  install.packages("corrplot")
library(corrplot)

if (! "precrec" %in% row.names(installed.packages()))
  install.packages("precrec")
library(precrec)


# ---------------------------------------------------------------------------------------------
# UTILS
# ---------------------------------------------------------------------------------------------
# default names of donor and acceptor tables
default_donors <- "donor.dat"
default_acceptors <- "acceptor.dat"

nucleotide_to_number <- function(nucleotide){
  if (nucleotide == "A"){
    return(1)
  }
  else if (nucleotide =="C"){
    return(2)
  }
  else if (nucleotide == "G"){
    return(3)
  }
  else if (nucleotide == "T"){
    return(4)
  }
  else if (nucleotide == "N")  {
    return(0)
  } else {
    warning("Not valid nucleotide -1 return as value.")
    return(-1)
  }
}

sequence_to_numbers <- function(sequence){
  return (map_dbl(unlist(strsplit(sequence,""), use.names=FALSE), nucleotide_to_number))
}

nucleotide_to_hoe <- function(nucleotide){
  if (nucleotide == "A"){
    return( c(1, 0, 0, 0))
  }
  else if (nucleotide =="C"){
    return( c(0, 1, 0, 0))
  }
  else if (nucleotide == "G"){
    return( c(0, 0, 1, 0))
  }
  else if (nucleotide == "T"){
    return( c(0, 0, 0, 1))
  }
  else if (nucleotide == "N")  {
    return( c(1, 1, 1, 1))
  } else {
    warning("Not valid nucleotide -1 return as value.")
    return(-1)
  }
}

sequence_to_hoe <- function(sequence){
  return (unlist(map(unlist(strsplit(sequence,""), use.names=FALSE), nucleotide_to_hoe)))
}


get_dataframe_from_sequences_dataframe_hoe <- function(file_name){
  data_table <- read.table(file_name)
  
  data_table <- matrix(data_table[,"V1"], nrow = 2, ncol=length(data_table[,"V1"])/2)
  data_table <- t(data_table)
  
  data_table <- cbind(data_table, do.call(rbind, map(data_table[, 2], sequence_to_hoe)))
  
  data_table <- cbind(data_table[, 1], data_table[, 3:ncol(data_table)])
  data_df <- as.data.frame(data_table)
  colnames(data_df)[1] <- c("class")
  remove(data_table)
  
  return(data_df)
}

get_dataframe_from_sequences_dataframe_num <- function(file_name){
  data_table <- read.table(file_name)
  
  data_table <- matrix(data_table[,"V1"], nrow = 2, ncol=length(data_table[,"V1"])/2)
  data_table <- t(data_table)
  
  
  data_table <- cbind(data_table, do.call(rbind, map(data_table[, 2], sequence_to_numbers)))
  
  data_table <- cbind(data_table[, 1], data_table[, 3:ncol(data_table)])
  data_df <- as.data.frame(data_table)
  colnames(data_df)[1] <- c("class")
  remove(data_table)
  
  return(data_df)
}



# ---------------------------------------------------------------------------------------------
# READ DATA
# ---------------------------------------------------------------------------------------------
acceptors_dataframe <- get_dataframe_from_sequences_dataframe_num(default_acceptors)
donors_dataframe <- get_dataframe_from_sequences_dataframe_num(default_donors)

features <- setdiff(names(acceptors_dataframe), "class")


# ---------------------------------------------------------------------------------------------
# REASERCH
# ---------------------------------------------------------------------------------------------
acceptors_cor = cor(acceptors_dataframe)
corrplot(acceptors_cor)


# ---------------------------------------------------------------------------------------------
# CLASSIFICATION
# ---------------------------------------------------------------------------------------------



# -----------------------
# DATA SPLIT
# -----------------------

split <- initial_split(acceptors_dataframe, prop = .5, strata = "class")
train_set <- training(split)
test_set  <- testing(split)

table(train_set$class) %>% prop.table()
table(test_set$class) %>% prop.table()

test_features <- test_set[, features]
test_class <- test_set$class

nb <- naiveBayes(class ~ + ., data=train_set)
pred <- predict(nb, test_features, type="raw")

precrec_obj <- evalmod(scores = pred[, 2], labels = test_class)
png(file = "plot.png", width=2000, height=1500, res=300)
autoplot(precrec_obj)
dev.off()
auc = auc(precrec_obj)
fileConn<-file("output.txt")
write.table(auc, fileConn)
close(fileConn)




tab <- table(pred, test_class)
caret::confusionMatrix(tab)

rf <- randomForest::randomForest(class ~ + ., data=train_set)
pred <- predict(rf, test_features)
tab <- table(pred, test_class)
caret::confusionMatrix(tab) 



# -----------------------
# CROSS VALIDATION
# -----------------------
x <- acceptors_dataframe[, features]
y <- acceptors_dataframe$class

train_control <- trainControl(
  method = "cv", 
  number = 10
)

# train model
nb_m1 <- train(
  x = x,
  y = y,
  method = "nb",
  trControl = train_control
)

# results
confusionMatrix(nb_m1)





