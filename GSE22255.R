data <- read.table("~/workspace/GSE22255/GSE22255_MAS5.MAS5.TXT", sep="\t", header = T)
lable <- names(data)
header.names <- data[,1]
data <- data[,-1]
data <- t(data)
data <- as.data.frame(data)