library(e1071)
load("~/workspace/GSE22255/GSE_DATA.RData")

cv <- function(data,feature){
  error.sum <- 0
  sets <- split(1:40,1:10)
  for (test in sets){
    
    training <- setdiff(1:nrow(data), test)
    
    model <- svm(data[training, feature], data$class[training] )
   
    pr <- predict(model, data[test,feature])
    #pr[pr > 0] <- 1
    #pr[pr < 0] <- 0
    cat("test set:", test, "true value:", convert.to.numeric(data$class[test]),"predicted value:" ,convert.to.numeric(pr), "\n")
    error.sum <- error.sum + length(which(!(pr == data$class[test]))) / 4
  }
  return(error.sum/10)
}

num.features <- ncol(data) - 1
error.rate <- 1.0
n <- 1
while(n < 100){
  two.features <- sample(num.features,2)
  cv.error.rate <- cv(data, two.features)
  
  if(cv.error.rate < error.rate){
    print(cv.error.rate)
    error.rate <- cv.error.rate
    cat("feature is:",two.features,", error rate is:",error.rate,"\n")
  }
  n <- n+1
  
}