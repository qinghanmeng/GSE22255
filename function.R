library(e1071)
library(randomForest)
load("~/workspace/GSE22255/GSE_DATA.RData")

get.tests.random <- function(x, nfold){
  x.random <- sample(x, length(x))
  sets <- split(x.random, 1:nfold)
  return(sets)
}

#for GSE22255
get.tests.balance <- function(){
  lable.one <- split(sample(1:20,20),1:10)
  lable.two <- split(sample(21:40,20),1:10)
  test <- rbind(as.data.frame(lable.one), as.data.frame(lable.two))
  return(test)
}

cv.svm <- function(x, y,nfold){
  error <- 0
  num.data <- nrow(x)
  #sets <- get.tests(1:num.data, nfold)
  sets <- get.tests.balance()
  for (test in sets){
    training <- setdiff(1:num.data, test)
    model <- svm(x[training,], y[training,] )
    pr <- predict(model, x[test,])
    cat("test set:", test, "true value:", convert.to.numeric(y[test,]),"predicted value:" ,convert.to.numeric(pr), "\n")
    error <- error + length(which(!(pr == y[test,])))
  }
  return(error/nrow(data))
}


# x, y are data frame type.
cv.logistic <- function(x, y, nfold){
  error <- 0
  num.data <- nrow(x)
  #sets <- get.tests(1:num.data, nfold)
  sets <- get.tests.balance()
  for (test in sets){
    training <- setdiff(1:num.data, test)
    training.set <- as.data.frame(cbind(x[training,],y[training,]))
    names(training.set)[length(training.set)] <- "class"
    model <- glm(class ~ ., training.set, family = binomial())
    pr <- predict(model, x[test,], type = "response")
    pr[pr > 0.5] <- 1
    pr[pr <= 0.5] <- 0
    cat("test set:", test, "true value:", convert.to.numeric(y[test,]),"predicted value:" ,convert.to.numeric(pr), "\n")
    error <- error + length(which(!(pr == y[test,])))
  }
  return(error/nrow(data))
}

print(cv.logistic(data[c(34456,21861)],data["class"] , 10))







