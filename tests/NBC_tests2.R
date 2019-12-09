if (! "e1071" %in% row.names(installed.packages()))
  install.packages("e1071")
library(e1071)

if (! "klaR" %in% row.names(installed.packages()))
  install.packages("klaR")
library(klaR)

if (! "caret" %in% row.names(installed.packages()))
  install.packages("caret")
library(caret)

# work, but values are bit stupid
test_2_model = data.frame(Class=c(1, 1, 1,  0, 0, 0), Len=c(2, 5, 5, 5, 2, 2) )

x = test_2_model[, 2, drop = FALSE]
y = test_2_model$Class
y = as.factor(y)


model = train(x,y,'nb')
predict(model, test_2_model[, c(2), drop=FALSE], type="prob")

model <- naiveBayes(x=x, y=y, laplace=1)
predict(model, test_2_model[, c(2), drop=FALSE], type="raw")

model <- klaR::NaiveBayes(x=x, grouping=y)
predict(model, test_2_model[, c(2), drop=FALSE], type="raw")





