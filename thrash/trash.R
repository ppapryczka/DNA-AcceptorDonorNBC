
#if (! "Hmisc" %in% row.names(installed.packages()))
#  install.packages("Hmisc")
#library(Hmisc)

#if (! "corrplot" %in% row.names(installed.packages()))
#  install.packages("corrplot")
#library(corrplot)

acceptors_cor = cor(acceptors_dataframe)
corrplot(acceptors_cor)
# stringr - for some string functions
#if (! "stringr" %in% row.names(installed.packages()))
#  install.packages("stringr")
#library(stringr)


# purrr - map function
if (! "purrr" %in% row.names(installed.packages()))
  install.packages("purrr")
library(purrr)