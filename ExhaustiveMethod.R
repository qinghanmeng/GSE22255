# this program is to combine two attributes using exhaustive method.
library(e1071)
library(iterators)
library(foreach)
load("~/workspace/GSE22255/GSE_DATA.RData")

cv.svm <- function(feature){
  error <- 0
  #num.data <- 40
  #sets <- get.tests(1:num.data, nfold)
  #sets <- get.tests.balance()
  sets <- split(1:40,1:10)
  
  for (test in sets){
    training <- setdiff(1:40, test)
    model <- svm(data[training,feature], data$class[training])
    pr <- predict(model, data[test,feature])
    #cat("test set:", test, "true value:", convert.to.numeric(y[test,]),"predicted value:" ,convert.to.numeric(pr), "\n")
    error <- error + length(which(!(pr == data$class[test])))
    
  }
  
  return(error/nrow(data))
}

test.per <- function(){
  for(i in 1:1000){
    cv.svm(svm.index)
  }
}


foreach(icount(1000)) %dopar% {
    cv.svm(svm.index)
}

#Timing stopped at: 246.644 0.14 247.307 
