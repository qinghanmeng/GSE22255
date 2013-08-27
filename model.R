library(e1071)
library(randomForest)

load("~/workspace/GSE22255/GSE_DATA.RData")

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


cv <- function(data,feature){
  #feature <- feature[1:8]
  error.sum <- 0
  for (i in seq(from = 1, to = 20, by = 2)){
    test <- c(i:(i + 1), (i+20):(i+21))
    training <- setdiff(1:nrow(data), test)
    #model <- naiveBayes(data[training, feature], data$class[training])
    model <- svm(data[training, feature], data$class[training] )
    
    #model <- glm(class ~ ., family = binomial(), data = data[training, c(feature,54676)])
    
    #model <- randomForest(data[training, feature], data$class[training],ntree = 5000)
    pr <- predict(model, data[test,feature])
    
    #pr[pr > 0] <- 1
    #pr[pr < 0] <- 2
    
    cat("test set", test, "true valeu:", data$class[test],"predicted value:",pr, "\n")
    error.sum <- error.sum + length(which(!(pr == data$class[test]))) / 4
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

#all.gene.features <- list(gene.1, gene.2, gene.3, gene.4, gene.5, gene.6, gene.7, gene.8)
#for (i in 1:(length(all.gene.features)-1)){
#  for(j in (i+1):(length(all.gene.features))){
#sub.features <- c(all.gene.features[[i]], all.gene.features[[j]])
#    cv(data, sub.features)
#}
#}

for( feature in feature_sets){
cv(data,feature)
}

#feature selection: t_statistics
#10-fold error rate: 0.375
#algorithm:naivebayes(GDA)
#select 100 features from each kind of ranking method, then intersect them and get following result
#[1]  8481  8489 30692  7515 26598 52341 22220 30636

#





