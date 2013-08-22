library(e1071)
library(randomForest)

#get features from file
info.gain.index <- read.table("~/workspace/rankgene-1.1/data/feature_selection/Information_Gain.list", sep = "\t", header = F)
info.gain.index <- info.gain.index$V1

gini.index <- read.table("~/workspace/rankgene-1.1/data/feature_selection/Gini_index.list", sep = "\t", header = F)
gini.index <- gini.index$V1

max.minority.index <- read.table("~/workspace/rankgene-1.1/data/feature_selection/Max_minority.list", sep = "\t", header = F)
max.minority.index <- max.minority.index$V1

sum.minority.index <- read.table("~/workspace/rankgene-1.1/data/feature_selection/Sum_minority.list", sep = "\t", header = F)
sum.minority.index <- sum.minority.index$V1

sum.variances.index <- read.table("~/workspace/rankgene-1.1/data/feature_selection/Sum_of_variances.list", sep = "\t", header = F)
sum.variances.index <- sum.variances.index$V1

t.statistic.index <- read.table("~/workspace/rankgene-1.1/data/feature_selection/t_statistic.list", sep = "\t", header = F)
t.statistic.index <- t.statistic.index$V1

twoing.rule.index <- read.table("~/workspace/rankgene-1.1/data/feature_selection/Twoing_rule.list", sep = "\t", header = F)
twoing.rule.index <- twoing.rule.index$V1

svm.index <- read.table("~/workspace/rankgene-1.1/data/feature_selection/svm.list", sep = "\t", header = F)
svm.index <- svm.index$V1


cv <- function(feature){
  feature <- feature[1:20]
  error.sum <- 0
  for (i in seq(from = 1, to = 20, by = 2)){
    test <- c(i:(i + 1), (i+20):(i+21))
    training <- setdiff(1:nrow(data), test)
    #model <- naiveBayes(data[training, feature], data$class[training])
    model <- svm(data[training, feature], data$class[training] )
    #model <- randomForest(data[training, feature], data$class[training])
    pr <- predict(model, data[test,feature])
    #cat("test set", test, "true valeu:", data$class[test],"predicted value:",pr, "\n")
    
    
    #print contingency table
    #print(table(pr,data$class[test]))
    error.sum <- error.sum + length(which(!(pr == data$class[test]))) / 4
  }
  print(error.sum/10)
}

feature_sets <- list(info.gain.index = info.gain.index,gini.index = gini.index,max.minority.index
                     = max.minority.index,sum.minority.index = sum.minority.index,
                     sum.variances.index = sum.variances.index,t.statistic.index = t.statistic.index,
                     twoing.rule.index = twoing.rule.index, svm.index = svm.index)

for( feature in feature_sets){
  cv(feature)
}


#feature selection: t_statistics
#10-fold error rate: 0.425
#algorithm:svm

#feature selection: t_statistics
#10-fold error rate: 0.375
#algorithm:naivebayes(GDA)
#select 100 features from each kind of ranking method, then intersect them and get following result
#[1]  8481  8489 30692  7515 26598 52341 22220 30636


