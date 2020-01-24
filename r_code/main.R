# modules - to split code
if (! "modules" %in% row.names(installed.packages()))
  install.packages("modules")
library(modules)

lib <- modules::use("R")