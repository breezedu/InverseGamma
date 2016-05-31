#read table
setwd("/Users/shuaiqizhang/Desktop/project /data ")
table<-read.table("sample1.txt")


#for the use of counting number of gene
sumenvarp<-aggregate(table$envarp, by=list(Category=table$gene), FUN=sum)
sumenvarpfc<-aggregate(table$envarpfc, by=list(Category=table$gene), FUN=sum)[,2]
tablesum<-data.frame(cbind(sumenvarp,sumenvarpfc))
colnames(tablesum)<-c("gene","sumenvarp","sumenvarpfc")