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


#construct random forest model and do cross validation
rf.cv <- function(feature){
  error.sum <- 0
  for (i in seq(from = 1, to = 20, by = 2)){
    test <- c(i:(i + 1), (i+20):(i+21))
    training <- setdiff(1:nrow(data), test)
    rf.model <- randomForest(data[training, feature], data$class[training])
    pr <- predict(rf.model, data[test,feature])
    error.sum <- error.sum + length(which(!(pr ==data$class[test]))) / 4
  }
  print(error.sum/10)
}

feature_sets <- list(info.gain.index = info.gain.index,gini.index = gini.index,max.minority.index
                     = max.minority.index,sum.minority.index = sum.minority.index,
                     sum.variances.index = sum.variances.index,t.statistic.index = t.statistic.index,
                     twoing.rule.index = twoing.rule.index, svm.index = svm.index)

for( feature in feature_sets){
  rf.cv(feature)
}