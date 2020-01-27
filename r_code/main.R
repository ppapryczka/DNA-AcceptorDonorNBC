# Główny skrypt do przeprowadzania eksperymentów.
# Autorzy: Patryk Pankiewicz, Łukasz Brzezicki

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)

source("utils.R")
source("donor_acceptors_NCB.R")

print(args)