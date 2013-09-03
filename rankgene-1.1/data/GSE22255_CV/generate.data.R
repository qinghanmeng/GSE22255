
test.sets <- split(1:40, 1:10)
lable <- cbind(as.character(raw.data[1,][3:42]), c(rep(0,20),rep(1,20)))

n <- 0
               
for(e in test.sets){
  training <- setdiff(1:42, e+2)
  write.table(lable[setdiff(1:40,e),], paste("~/workspace/GSE22255/rankgene-1.1/data/GSE22255_CV/",n,"/",n,".lable", sep = ""), sep = "\t",, quote = F, col.names = F, row.names = F)
  #write.table(raw.data[,training], paste("~/workspace/GSE22255/rankgene-1.1/data/GSE22255_CV/",n,"/",n, sep = ""), sep = "\t",, quote = F, col.names = F, row.names = F)
  n <- n + 1
}