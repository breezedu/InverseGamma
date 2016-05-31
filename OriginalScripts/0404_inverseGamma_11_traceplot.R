## 0404_inverseGamma_11.R

#setwd("/Users/shuaiqizhang/Desktop/project /data ")
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


#############################################################################################
#############################################################################################


hiernormalinvg<-"
data{ #get the data set 
	int<lower=0> N;   		# number of exon level
	int<lower=0> J;   		# number of gene level
	int <lower=1,upper=J> gene[N];  #index of gene
	int <lower=1,upper=N> exon[N];  #index of exon 
	vector[N] xij;   		#x at exon level
	vector[N] yij; 			#y at exon level
}
parameters{ 				#specify the parameter we want to know 
	real beta;  			#common slope for the exon level
	real mu;      			#common intercept for the exon level
	vector[N] aij; 			#random intercept for the exon level
	real <lower=0> sigma_aj2[J];  	#variance of intercept at exon level 
	vector[J] aj; 			#random intercept for the gene level 
	real <lower=0> sigma_a;  	#variance of intercept at gene level
	real <lower=0> sigma; 		#variance of yij
}
transformed parameters{ #specify the model we will use 
	
}
model { #give the prior distribution 
   vector[N] lambdaij; #exon level model
    for (i in 1:N)
        lambdaij[i] <- mu+beta*xij[i]+aij[exon[i]]+aj[gene[i]];#specify the gene group
   beta ~normal(0,50);
   mu~normal(0,50);
   sigma ~ uniform(0, 20);
   sigma_a ~uniform(0,10);
   sigma_aj2 ~inv_gamma(1,1);
   aj ~ normal(0, sigma_a);
   yij ~ normal(lambdaij,sigma);
   for (i in 1:N)
       aij[i]~ normal(0,sqrt(sigma_aj2[gene[i]]));
   }

"


library("rstan")

J<-dim(table1)[1] #gene number 
N<-dim(table)[1]  #exon number 
xij=c(table$envarp)
yij=log(table$envarpfc+0.01)

gene<-as.numeric(table$gene) #list
genelevel<-length(unique(gene)) #number
indexg<-match(gene, unique(gene))  #list
exon<-c(1:length(table$envarpfc))
M1_table<-list(N=N, J=J, xij=xij,
yij=yij,gene=indexg, exon=exon)
control=list(adapt_delta=0.99,max_treedepth=12)
fitinv<-stan(model_code=hiernormalinvg, data=M1_table,iter=40000,warmup=35000,chains=4) #10,000 samplings 
#fitinv<-stan(model_code=hiernormalinvg, data=M1_table,iter=100,chains=4) #10,000 samplings 
answer<-extract(fitinv,permuted=TRUE)

pdf(file = "0404_inverseGamma_11_sampleTracePlot.pdf")
traceplot(fitinv,par="mu")
traceplot(fitinv,par="beta")
traceplot(fitinv,par="sigma")

print(fitinv)


#########################################################
## END
#########################################################