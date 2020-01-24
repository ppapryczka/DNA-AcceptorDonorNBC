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

# ---------------------------------------------------------------------------------------------
# CREATE FUNCTIONS AND DECLARE CONSTS
# ---------------------------------------------------------------------------------------------

# default names of donor and acceptor tables
default_donors <- "donor.dat"
default_acceptors <- "acceptor.dat"

# function that converts nucleotide to number
#nucleotide_to_number <- function(nucleotide){
#  if (nucleotide == "A"){
#    return( c(1, 0, 0, 0))
#  }
#  else if (nucleotide =="C"){
#    return( c(0, 1, 0, 0))
#  }
#  else if (nucleotide == "G"){
#    return( c(0, 0, 1, 0))
#  }
#  else if (nucleotide == "T"){
#    return( c(0, 0, 0, 1))
#  }
#  else if (nucleotide == "N")  {
#    return( c(1, 1, 1, 1))
#  } else {
#    warning("Not valid nucleotide -1 return as value.")
#    return(-1)
#  }
#}

# function that converts nucleotide to number
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

#sequence_to_numbers <- function(sequence){
#  return (unlist(map(unlist(strsplit(sequence,""), use.names=FALSE), nucleotide_to_number)))
#}

get_columns_name_for_DNA_one_hot_encoding <- function(length){
  l = c("class")
  for(n in c(0:length)){
    l <- c(l, paste(toString(n), "_a", sep=""), paste(toString(n), "_c", sep=""), paste(toString(n), "_g", sep=""), paste(toString(n), "_t", sep=""))
    
  }
  return (l)
}

# ---------------------------------------------------------------------------------------------
# READ DATA
# ---------------------------------------------------------------------------------------------

# read data from tables
donors_table <- read.table(default_donors)
acceptors_table <- read.table(default_acceptors)


# convert acceptors table
acceptors_table <- matrix(acceptors_table[,"V1"], nrow = 2, ncol=length(acceptors_table[,"V1"])/2)
acceptors_table <- t(acceptors_table)

# convert donors table
donors_table <- matrix(donors_table[,"V1"], nrow = 2, ncol=length(donors_table[,"V1"])/2)
donors_table <- t(donors_table)

# get columns names for dataframe
seq_length = nchar(acceptors_table[1, 2])
df_columns = get_columns_name_for_DNA_one_hot_encoding(seq_length-1)

# ---------------------------------------------------------------------------------------------
# PRODUCE ATTRIBUTES
# ---------------------------------------------------------------------------------------------

# explonation of calling <- map_dbl <- sequence_to_numbers <- map <- do.call/rbind <- cbind

# add number-position arguments to acceptors table
acceptors_table <- cbind(acceptors_table, do.call(rbind, map(acceptors_table[, 2], sequence_to_numbers)))

# add number-position arguments to donors table
donors_table <- cbind(donors_table, do.call(rbind, map(donors_table[, 2], sequence_to_numbers)))

# get rid of "seq" fields
acceptors_table = cbind(acceptors_table[, 1], acceptors_table[, 3:ncol(acceptors_table)]) 
donors_table = cbind(donors_table[,1 ], donors_table[, 3:ncol(donors_table)])


# ---------------------------------------------------------------------------------------------
# ANALYSIS AND MODEL
# ---------------------------------------------------------------------------------------------

# convert both tables to dataframes and delete tables
acceptors_dataframe = as.data.frame(acceptors_table)
colnames(acceptors_dataframe)[1] <- c("class")
#colnames(acceptors_dataframe) <- df_columns
remove(acceptors_table)

donors_dataframe = as.data.frame(donors_table)
colnames(donors_dataframe)[1] <- c("class")
#colnames(donors_dataframe) <- df_columns
remove(donors_table)


#summary(acceptors_dataframe[acceptors_dataframe[, 1]==1, ,drop=FALSE])
#summary(acceptors_dataframe[acceptors_dataframe[, 1]==0, ,drop=FALSE])
#summary(acceptors_dataframe[acceptors_dataframe[, 1]==0, ,drop=FALSE][, c(1, 5)])

# ---------------------------------------------------------------------------------------------
# CLASSIFICATION
# ---------------------------------------------------------------------------------------------

#model <- naiveBayes(Class ~ ., data = acceptors_dataframe[, c(1:10),])
#predict(model, acceptors_dataframe[,3: ncol(acceptors_dataframe)])
#predict(model, acceptors_dataframe[,3: ncol(acceptors_dataframe)], type = "raw")

#pred <- predict(model, acceptors_dataframe)
#table(pred, acceptors_dataframe$Class)


#x <- acceptors_dataframe[, , drop=FALSE]
#y <- acceptors_dataframe$Class
#y <- as.factor(y)
#model <- naiveBayes(x=x, y=y, laplace = 0)
#pred <- predict(model, acceptors_dataframe[,c(2,3,4,5), drop=FALSE])
#predict(model, acceptors_dataframe[,c(2,3,4,5), drop=FALSE], type = "raw")
#table(pred, acceptors_dataframe$Class)

split_index = sort(sample(nrow(acceptors_dataframe), nrow(acceptors_dataframe)*0.5))
train <- acceptors_dataframe[split_index,]
test <- acceptors_dataframe[-split_index,]
features <- setdiff(names(acceptors_dataframe), "class")

y_test = test$class
test = test[, features]

#features <- features[sample(ncol(acceptors_dataframe) - 1, 3)]
#y_test = test$class
#test = test[, features]
#train = train[, c("class", features)]

nb <- naiveBayes(class ~ + ., data=train)
pred <- predict(nb, test)
tab <- table(pred, y_test)
caret::confusionMatrix(tab)  

rf <- randomForest::randomForest(class ~ + ., data=train)
pred <- predict(rf, test)
tab <- table(pred, y_test)
caret::confusionMatrix(tab)  

data("HouseVotes84")


model <- naiveBayes(Class ~., data=HouseVotes84)
pred <- predict(model, HouseVotes84[1: 10,])
table(pred, HouseVotes84[1:10,]$Class)

split <- initial_split(acceptors_dataframe, prop = .5, strata = "class")
train <- training(split)
test  <- testing(split)

table(train$class) %>% prop.table()
table(test$class) %>% prop.table()

features <- setdiff(names(acceptors_dataframe), "class")
x <- acceptors_dataframe[, features]
y <- acceptors_dataframe$class

train_control <- trainControl(
  method = "cv", 
  number = 10
)

nb.m1 <- train(
  x = x,
  y = y,
  method = "nb",
  trControl = train_control
)

confusionMatrix(nb.m1)

data("iris")

nb <- naiveBayes(Species ~., data=iris)
pred <- predict(nb, iris) 
tab <- table(pred, iris$Species)
caret::confusionMatrix(tab)  

split <- initial_split(iris, prop = .7, strata = "Species")
train <- training(split)
test  <- testing(split)

table(train$Species) %>% prop.table()
table(test$Species) %>% prop.table()

features <- setdiff(names(train), "Species")
x <- train[, features]
y <- train$Species

train_control <- trainControl(
  method = "cv", 
  number = 10
)
search_grid <- expand.grid(
  usekernel = c(TRUE, FALSE),
  fL = 0:5,
  adjust = seq(0, 5, by = 1)
)


nb.m1 <- train(
  x = x,
  y = y,
  method = "nb",
  tune_grid= search_grid,
  trControl = train_control
)

confusionMatrix(nb.m1)


# ---------------------------------------------------------------------------------------------
# RESULTS
# ---------------------------------------------------------------------------------------------
