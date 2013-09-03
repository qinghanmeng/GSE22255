library(e1071)
library(randomForest)


#get features from file

info.gain.index <- read.table("~/workspace/GSE22255/rankgene-1.1/data/feature_selection/Information_Gain.list", sep = "\t", header = F)
info.gain.index <- info.gain.index$V1

gini.index <- read.table("~/workspace/GSE22255/rankgene-1.1/data/feature_selection/Gini_index.list", sep = "\t", header = F)
gini.index <- gini.index$V1

max.minority.index <- read.table("~/workspace/GSE22255/rankgene-1.1/data/feature_selection/Max_minority.list", sep = "\t", header = F)
max.minority.index <- max.minority.index$V1

sum.minority.index <- read.table("~/workspace/GSE22255/rankgene-1.1/data/feature_selection/Sum_minority.list", sep = "\t", header = F)
sum.minority.index <- sum.minority.index$V1

sum.variances.index <- read.table("~/workspace/GSE22255/rankgene-1.1/data/feature_selection/Sum_of_variances.list", sep = "\t", header = F)
sum.variances.index <- sum.variances.index$V1

t.statistic.index <- read.table("~/workspace/GSE22255/rankgene-1.1/data/feature_selection/t_statistic.list", sep = "\t", header = F)
t.statistic.index <- t.statistic.index$V1

twoing.rule.index <- read.table("~/workspace/GSE22255/rankgene-1.1/data/feature_selection/Twoing_rule.list", sep = "\t", header = F)
twoing.rule.index <- twoing.rule.index$V1

svm.index <- read.table("~/workspace/GSE22255/rankgene-1.1/data/feature_selection/svm.list", sep = "\t", header = F)
svm.index <- svm.index$V1

convert.to.numeric <-function(f){
  return(as.numeric(as.character(f)))
}

cv <- function(data,feature,m.name){
  sets <- split(1:nrow(data),1:9)
  error.sum <- 0
  for (test in sets){
    
    training <- setdiff(1:nrow(data), test)
    
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
    error.sum <- error.sum + length(which(!(pr == data$class[test]))) / 7
  }
  print(error.sum/10)
}



feature_sets <- list(info.gain.index = info.gain.index,gini.index = gini.index,max.minority.index
                     = max.minority.index,sum.minority.index = sum.minority.index,
                     sum.variances.index = sum.variances.index,t.statistic.index = t.statistic.index,
                     twoing.rule.index = twoing.rule.index, svm.index = svm.index)


normalize.data <- function(data){
  data.normal <- lapply(data[,-ncol(data)], function(e) (e-mean(e) /sd(e)))
  data.normal$class <- data$class
  return(data.normal)
}

#data.normal <- as.data.frame(normalize.data(data))

all.gene.features <- list(gene.1, gene.2, gene.3, gene.4, gene.5, gene.6, gene.7, gene.8)

combine.gene <- function(){
  for (i in 1:(length(all.gene.features)-1)){
    for(j in (i+1):(length(all.gene.features))){
      sub.features <- c(all.gene.features[[i]], all.gene.features[[j]])
      cv(data, sub.features)
    }
  }
}

#for( feature in feature_sets){
#  cv(GSE16561,GSE16561.svm.feature,"svm")
#}
cv(GSE16561,GSE16561.svm.feature[1:10],"svm")





