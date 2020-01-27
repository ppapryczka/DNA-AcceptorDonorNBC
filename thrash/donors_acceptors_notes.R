# old file, now just notes, maybe use some of attributes  

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
columns_ids = c()

for ( x in 1:3 )
{
  prm <-permutations(n = 4, r = x,  c("A", "G", "T", "C"), repeats.allowed = TRUE)
  prm <- apply(prm, 1, function(x)paste0(x, collapse=''))
  columns_ids <- c(columns_ids, prm)
  
  # count occurances of each permutation in sequences
  for ( p in prm )
  {
    acceptors_table <- cbind(acceptors_table, matrix(str_count(acceptors_table[, 2], p)))
    donors_table <- cbind(donors_table, matrix(str_count(donors_table[, 2], p)))
  }
}

for ( x in 1:3)
{
  prm <-permutations(n = 5, r = x,  c("A", "G", "T", "C", "."))
  prm <- apply(prm, 1, function(x)paste0(x, collapse=''))
  columns_ids <- c(columns_ids, prm)
  
  #count occurances of each permutation in sequences
  for ( p in prm )
  {
    acceptors_table <- cbind(acceptors_table, matrix(length(unlist(str_extract_all(acceptors_table[, 2], p))))) 
    donors_table <- cbind(donors_table, matrix(length(unlist(str_extract_all(donors_table[, 2], p)))))
  }
}

for (x in acceptors_table[, 2])
{
  print(unlist(str_extract_all(x, ".")))
  print("Length")
  print(length(unlist(str_extract_all(x, ".."))))
  print(str_count(x, "."))
  print(str_count(x, ".."))
}


summary(acceptors_table)
print(columns_ids)





# ---------------------------------------------------------------------------------------------
# DUMMY
# ---------------------------------------------------------------------------------------------

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
