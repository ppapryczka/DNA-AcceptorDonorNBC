if (! "stringr" %in% row.names(installed.packages()))
  install.packages("stringr")

library(stringr)

if (! "gtools" %in% row.names(installed.packages()))
  install.packages("gtools")

library(gtools)


# default names of donor and acceptor tables
default_donors <- "donor.dat"
default_acceptors <- "acceptor.dat"

# read data from tables
donors_table <- read.table(default_donors)
acceptors_table <- read.table(default_acceptors)


# convert acceptors table
acceptors_table <- matrix(acceptors_table[,"V1"], nrow = 2, ncol=length(acceptors_table[,"V1"])/2)
acceptors_table <- t(acceptors_table)

# convert donors table
donors_table <- matrix(donors_table[,"V1"], nrow = 2, ncol=length(donors_table[,"V1"])/2)
donors_table <- t(donors_table)

# get all permutations of A G T C nucleotides of specific length 
for ( x in 1:4 )
{
  prm <-permutations(n = 4, r = x,  c("A", "G", "T", "C"))
  prm <- apply(prm, 1, function(x)paste0(x, collapse=''))
  
  # count occurances of each permutation in sequences
  for ( p in prm )
  {
    acceptors_table <- cbind(acceptors_table, matrix(str_count(acceptors_table[, 2], p)))
    donors_table <- cbind(donors_table, matrix(str_count(donors_table[, 2], p)))
  }
}


