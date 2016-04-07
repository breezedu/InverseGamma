setwd("/Users/shuaiqizhang/Desktop/project /data ")
table<-read.table("exon_level_process_v2.txt")
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
table$x=scale(table$envarp)

table1<-table[8700:9200,]
table1<-table1[which(table1$envarp!=0), ]
#get rid of genes with only one exon
g<-unique(table1$gene)
gene<-as.numeric(table1$gene) #list
genelevel<-length(unique(gene)) #number

indexg<-sort(match(gene, unique(gene)) )#list
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

delet=which(indexg%in%which(leng==1)) #get the index of exon in gene which only contains one exon
table1=table1[-delet,]# delet the row of genes contains only one exons 
#table1[which(table2$gene=="GNB1"),]



table2<-table[139700:140200,]
table2<-table2[which(table2$envarp!=0), ]

#get rid of genes with only one exon
g<-unique(table2$gene)
gene<-as.numeric(table2$gene) #list
genelevel<-length(unique(gene)) #number

indexg<-sort(match(gene, unique(gene)) )#list
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
leng=c(leng,len)  #has problem for the whole data set careful

delet=which(indexg%in%which(leng==1)) #get the index of exon in gene which only contains one exon
table2=table2[-delet,]# delet the row of genes contains only one exons 
#table2[which(table2$gene=="ABCC9"),]

table<-rbind(table1,table2) #828 exons

#for the use of counting number of gene
sumenvarp<-aggregate(table$envarp, by=list(Category=table$gene), FUN=sum)
sumenvarpfc<-aggregate(table$envarpfc, by=list(Category=table$gene), FUN=sum)[,2]
tablesum<-data.frame(cbind(sumenvarp,sumenvarpfc))
colnames(tablesum)<-c("gene","sumenvarp","sumenvarpfc")


