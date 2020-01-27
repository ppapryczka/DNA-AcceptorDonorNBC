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
png(file = "plot.png", width=2000, height=1500, res=300)
autoplot(precrec_obj)
dev.off()
auc = auc(precrec_obj)

file.create("dummy.txt")
write.table(auc, "dummy.txt", append=TRUE)

pred <-predict(nb, test_features)
tab <- table(pred, test_class)
caret::confusionMatrix(tab)
write.table(as.table(confusionMatrix(tab)), "dummy.txt", append=TRUE)


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





