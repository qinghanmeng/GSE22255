library(e1071)
library(randomForest)

for(i in 0:9){
  assign(paste("feature.", i+1, sep = ""),read.table(paste("~/workspace/GSE22255/rankgene-1.1/data/GSE22255_CV/",i,"/",i,".index",sep="")))
  temp <- get(paste("feature.", i+1, sep = ""))
  assign(paste("feature.", i+1, sep = ""), temp$V1)
}
convert.to.numeric <-function(f){
  return(as.numeric(as.character(f)))
}


features <- list()

for(i in 1:10){
  features <- cbind(features, get(paste("feature.",i,sep = "")))
}

cv <- function(data,feature.sets,m.name){
  error.sum <- 0
  sets <- split(1:40,1:10)
  for(i in 1:10){
    test <- sets[[i]]
    training <- setdiff(1:40, test)
    
    feature <- as.integer(feature.sets[,i])
  
    model <- NULL
    if(m.name == "naiveBayes"){
      model <- naiveBayes(data[training, feature], data$class[training])
    }
    
    if(m.name == "svm"){
      model <- svm(data[training, feature], data$class[training] )
    }
    
    if(m.name == "logistic"){
      model <- glm(class ~ ., family = binomial(), data = data[training, c(feature,54676)])
    }
    
    if(m.name == "randomForest"){
      model <- randomForest(data[training, feature], data$class[training],ntree = 5000)
    }
    
    pr <- predict(model, data[test,feature])
    #pr[pr > 0] <- 1
    #pr[pr < 0] <- 2
    cat("test set:", test, "true value:", convert.to.numeric(data$class[test]),"predicted value:" ,convert.to.numeric(pr), "\n")
    error.sum <- error.sum + length(which(!(pr == data$class[test]))) / 4
  }
  print(error.sum/10)
}

cv(data,features,"svm")
