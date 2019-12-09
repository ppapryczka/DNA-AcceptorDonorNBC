# ---------------------------------------------------------------------------------------------
# LOAD PACKAGES 
# ---------------------------------------------------------------------------------------------

# stringr - for some string functions
if (! "stringr" %in% row.names(installed.packages()))
  install.packages("stringr")
library(stringr)

# gtools - permutation and combinations functions
if (! "gtools" %in% row.names(installed.packages()))
  install.packages("gtools")
library(gtools)

# purrr - map function
if (! "purrr" %in% row.names(installed.packages()))
  install.packages("purrr")
library(purrr)

# e1071 - naive bayes
if (! "e1071" %in% row.names(installed.packages()))
  install.packages("e1071")
library(e1071)

# ---------------------------------------------------------------------------------------------
# CREATE FUNCTIONS AND DECLARE CONSTS
# ---------------------------------------------------------------------------------------------

# default names of donor and acceptor tables
default_donors <- "donor.dat"
default_acceptors <- "acceptor.dat"

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
colnames(acceptors_dataframe)[1] <- c("Class")
remove(acceptors_table)

donors_dataframe = as.data.frame(donors_table)
colnames(donors_dataframe)[1] <- c("Class")
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


x <- acceptors_dataframe[, c(2, 3, 4, 5), drop=FALSE]
y <- acceptors_dataframe$Class
y <- as.factor(y)
model <- naiveBayes(x=x, y=y, laplace = 0)
pred <- predict(model, acceptors_dataframe[,c(2,3,4,5), drop=FALSE])
predict(model, acceptors_dataframe[,c(2,3,4,5), drop=FALSE], type = "raw")
table(pred, acceptors_dataframe$Class)


# ---------------------------------------------------------------------------------------------
# RESULTS
# ---------------------------------------------------------------------------------------------
