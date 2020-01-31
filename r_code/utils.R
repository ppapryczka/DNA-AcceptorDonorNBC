# Narzędzia do wczytywania i przetwarzania sekwencji DNA na atrybuty.
# Autorzy: Patryk Pankiewicz, Łukasz Brzezicki

# ---------------------------------------------------------------------------------------------
# LOAD PACKAGES 
# ---------------------------------------------------------------------------------------------

# purrr - map functions
if (! "purrr" %in% row.names(installed.packages()))
  install.packages("purrr")
library(purrr)


# ---------------------------------------------------------------------------------------------
# FUNCTIONS 
# ---------------------------------------------------------------------------------------------

#' Convert given nucleotide to number A=1, C=2, G=3, T=4.
#' @param nucleotide Nucleotide to convert
#' @return Code of nucleotide.
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

#' Convert sequence of nucleotides to matix of numbers
#' @param sequence Sequence of nucleotides to convert
#' @return Matrix of codes of nucelotides.
sequence_to_numbers <- function(sequence){
  return (map_dbl(unlist(strsplit(toString(sequence),""), use.names=FALSE), nucleotide_to_number))
}

#' Convert nucleotide to hot-one-enconding.
#' @param nucleotide Nucleotide to convert
#' @return Hot-one-encoding list got nucelotidee.
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

#' Convert sequence of nucleotides to matix of numbers
#' @param sequence Sequence of nucleotides to convert
#' @return Matrix of codes of nucelotides.
sequence_to_hoe <- function(sequence){
  return (unlist(map(unlist(strsplit(sequence,""), use.names=FALSE), nucleotide_to_hoe)))
}

#' Load data from given file and sonvert sequences from file using one-hot-encoding.
#' @param file_name Path to file name.
#' @return Dataframe with sequences coded in one_hot_encoding.
get_dataframe_from_sequences_dataframe_hoe <- function(file_name){
  data_table <- read.table(file_name)
  
  data_table <- matrix(data_table[,"V1"], nrow = 2, ncol=length(data_table[,"V1"])/2)
  data_table <- t(data_table)
  
  data_table <- cbind(data_table, do.call(rbind, map(data_table[, 2], sequence_to_hoe)))
  
  data_table = cbind(data_table[, 1], data_table[, 3:ncol(data_table)])
  data_df = as.data.frame(data_table)
  colnames(data_df)[1] <- c("class")
  remove(data_table)
  
  return(data_df)
}

#' Load data from given file and sonvert sequences from file using one-hot-encoding.
#' @param file_name Path to file name.
#' @return Dataframe with sequences coded in one_hot_encoding.
get_dataframe_from_sequences_dataframe_num <- function(file_name){
  data_table <- read.table(file_name)
  
  data_table <- matrix(data_table[,"V1"], nrow = 2, ncol=length(data_table[,"V1"])/2)
  data_table <- t(data_table)
  
  data_table <- cbind(data_table, do.call(rbind, map(data_table[, 2], sequence_to_numbers)))
  
  data_table = cbind(data_table[, 1], data_table[, 3:ncol(data_table)])
  data_df = as.data.frame(data_table)
  colnames(data_df)[1] <- c("class")
  remove(data_table)
  
  return(data_df)
}
