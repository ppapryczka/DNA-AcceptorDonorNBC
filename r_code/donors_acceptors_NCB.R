#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# ---------------------------------------------------------------------------------------------
# LOAD PACKAGES 
# ---------------------------------------------------------------------------------------------

# e1071 - naive bayes
if (! "e1071" %in% row.names(installed.packages()))
  install.packages("e1071")
library(e1071)

# caret - confusion matrix
if (! "caret" %in% row.names(installed.packages()))
  install.packages("caret")
library(caret)

# randomForest - ... random forest
if (! "randomForest" %in% row.names(installed.packages()))
  install.packages("randomForest")
library(randomForest)

# datarsample - data splitig
if (! "rsample" %in% row.names(installed.packages()))
  install.packages("rsample")
library(rsample)

# precrec - pr and auc curve 
if (! "precrec" %in% row.names(installed.packages()))
  install.packages("precrec")
library(precrec)


default_donors <- "donor.dat"
default_acceptors <- "acceptor.dat"

source("utils.R")


# ---------------------------------------------------------------------------------------------
# FUNCTIONS 
# ---------------------------------------------------------------------------------------------

create_report_file <- function(pred_class, pred_raw, true_class, report_path, prefix =""){
  # save roc curve and pr curve
  plot_path = file.path(report_path, paste(prefix, "plot.png", sep = "_"))
  precrec_obj <- evalmod(scores = pred_raw[, 2], labels = true_class)
  png(file = plot_path, width=3000, height=1500, res=300)
  autoplot(precrec_obj)
  dev.off()
  
  # create raport file
  report_file_path = file.path(report_path, paste(prefix, "report.txt", sep="_"))
  file.create(report_file_path)
  
  # save aucs
  sink(report_file_path)
  auc = auc(precrec_obj)
  print(auc)
  
  # save confusion matrix
  tab <- table(pred_class, true_class)
  print(caret::confusionMatrix(tab))
  sink()
}

train_test_validation_nb <- function(df, split_prop = 0.5, report_path ="."){
  # get only features names
  features <- setdiff(names(df), "class")
  
  # get split 
  split <- initial_split(df, prop = split_prop, strata = "class")
  
  # split data
  train_set <- training(split)
  test_set  <- testing(split)
  
  # split features and target class
  test_features <- test_set[, features]
  test_class <- test_set$class
  
  # train classifier
  nb <- naiveBayes(class ~ + ., data=train_set)
  
  # predict raw 
  pred_raw <- predict(nb, test_features, type="raw")
  # predict class
  pred_class <-predict(nb, test_features)
  
  create_report_file(pred_class, pred_raw, test_class, report_path, "nb")
}

train_test_validation_rf <- function(df, split_prop = 0.5, report_path ="."){
  # get only features names
  features <- setdiff(names(df), "class")
  
  # get split 
  split <- initial_split(df, prop = split_prop, strata = "class")
  
  # split data
  train_set <- training(split)
  test_set  <- testing(split)
  
  # split features and target class
  test_features <- test_set[, features]
  test_class <- test_set$class
  
  # train classifier
  rf <- randomForest(class ~ + ., data=train_set)
  
  # predict raw 
  pred_raw <- predict(rf, test_features, type="prob")
  # predict class
  pred_class <-predict(rf, test_features)
  
  create_report_file(pred_class, pred_raw, test_class, report_path, "rf")
}

cross_validation <- function(df, report_path =".", method = "nb", cv_num = 10){
  features <- setdiff(names(df), "class")
  
  # split class and features
  x <- df[, features]
  y <- df$class
  
  # init train control 
  train_control <- trainControl(
    method = "cv", 
    number = cv_num
  )
  
  # train model
  clf <- train(
    x = x,
    y = y,
    method = method,
    trControl = train_control
  )

  # predict raw  
  pred_raw <- predict(clf, type="prob")
  # predict class
  pred_class <- predict(clf, type="raw")
  
  create_report_file(pred_class, pred_raw, y, report_path, method)
}

df = get_dataframe_from_sequences_dataframe_num(default_acceptors)
train_test_validation_nb(df, report_path = "./dupa")

df = get_dataframe_from_sequences_dataframe_hoe(default_acceptors)
train_test_validation_nb(df, report_path = ".")

df = get_dataframe_from_sequences_dataframe_num(default_acceptors)
cross_validation(df, report_path = "cv", method = "nb")



# ---------------------------------------------------------------------------------------------
# READ DATA
# ---------------------------------------------------------------------------------------------
acceptors_dataframe <- get_dataframe_from_sequences_dataframe_num(default_acceptors)
donors_dataframe <- get_dataframe_from_sequences_dataframe_num(default_donors)

features <- setdiff(names(acceptors_dataframe), "class")


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
png(file = "plot.png", width=3000, height=1500, res=300)
autoplot(precrec_obj)
dev.off()
auc = auc(precrec_obj)
print(auc)


pred <-predict(nb, test_features)
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
train()

# results
confusionMatrix(nb_m1)
pred_prob <- predict(nb_m1, type="prob")
pred_raw <- predict(nb_m1, type="raw")
tab <- table(pred_raw, y)
caret::confusionMatrix(tab)



features <- setdiff(names(df), "class")

# split class and features
x <- acceptors_dataframe[, features]
y <- acceptors_dataframe$class

# init train control 
train_control <- trainControl(
  method = "cv", 
  number = 10
)

# train model
clf <- train(
  x = x,
  y = y,
  method = "nb",
  trControl = train_control
)

# predict raw  
pred_raw <- predict(clf, type="prob")
# predict class
pred_class <- predict(clf, type="raw")


print(caret::confusionMatrix(tab))

precrec_obj <- evalmod(scores = pred_raw[, 2], labels = y)
autoplot(precrec_obj)

