

setwd("/Users/shuaiqizhang/Desktop/project /data ")
table<-read.table("/dscrhome/gd44/SQProject/RStan/2016_spring/exon_level_process_v2.txt")
colnames(table)<- c("chr", "gene", "dom", "subdom", "exon", "gene.dom", 
             "gene.dom.subdom", 
             "envarp",    # pass
             "envarpf",   # pass functional
             "envarpfr",  # pass functional rare
             "emutr")     # mutation rate 

table<-within(table,envarpfc<-envarpf-envarpfr)#y
table<-within(table,gene<-factor(gene))
table<-within(table,gene.dom<-factor(gene.dom))
table<-within(table,gene.dom.subdom<-factor(gene.dom.subdom))
table<-table[which(table$envarp!=0), ]
table$x=scale(table$envarp)
table<-table[1:1000,]

#get rid of genes with only one exon
gene<-as.numeric(table$gene) #list
genelevel<-length(unique(gene)) #number
indexg<-match(gene, unique(gene))  #list
leng=c()
len=1
for  ( i in 1:(length(indexg)-1)){
	if (indexg[i+1]==indexg[i]){
		len=len+1
	   }
	if (indexg[i+1]!=indexg[i]){
		leng=c(leng,len)
		len=1
	}
}
leng=c(leng,len)   #has problem for the whole data set careful

delet=which(indexg%in%which(leng%in%1)) #get the index of exon in gene which only contains one exon
table=table[-delet,]# delet the row of genes contains only one exons 


#for the use of counting number of gene
sumenvarp<-aggregate(table$envarp, by=list(Category=table$gene), FUN=sum)
sumenvarpfc<-aggregate(table$envarpfc, by=list(Category=table$gene), FUN=sum)[,2]
table1<-data.frame(cbind(sumenvarp,sumenvarpfc))
colnames(table1)<-c("gene","sumenvarp","sumenvarpfc")


