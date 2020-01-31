#!/usr/bin/env Rscript

# Główny skrypt do przeprowadzania eksperymentów.
# Autorzy: Patryk Pankiewicz, Łukasz Brzezicki

args = commandArgs(trailingOnly=TRUE)

# check arguments
if (length(args) < 3){
  cat("Arguments error!\n")
  cat("Usage: <path_to_data_file> <path_to_report_folder> <validation method>\n")
  cat("Validation metods:  'sv' - simple validation, 'cv'- cross validation\n.")
  quit(save = "no", status = 0)
}

if (!identical(args[3], "sv") && !identical(args[3], "cv")){
  cat("Wrong validation method given!\n")
  cat("Validation metods:  'sv' - simple validation, 'cv'- cross validation\n.")
  quit(save = "no", status = 0)
}

# load other files
source("utils.R")
source("donors_acceptors_NCB.R")

# run classification

if(identical(args[3], "sv")){
  df = get_dataframe_from_sequences_dataframe_num(args[1])
  print(summary(df))
  train_test_validation_nb(df, report_path = args[2])  
}

if(identical(args[3], "cv")){
  df = get_dataframe_from_sequences_dataframe_num(args[1])
  print(summary(df))
  cross_validation(df, report_path = args[2])
}

quit(save = "no", status = 0)
  