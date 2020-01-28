args <- commandArgs(trailingOnly = TRUE)
x <- as.integer(args[1])
y <- as.integer(args[2])

cat(x, y)source("utils.R")
source("donors_acceptors_NCB.R")

format(Sys.time(), "%Y_%m_%d_%H_%M_%S")

strptime(Sys.time(),format='%Y-%m-%d %H:%M:%S')

quit(save = "no", status = 0)
print(strptime(Sys.time(), ))