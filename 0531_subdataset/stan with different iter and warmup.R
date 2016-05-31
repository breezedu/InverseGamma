
## 0531 RStan with subdataset

#read table
##setwd("/Users/shuaiqizhang/Desktop/project /data ")
table<-read.table("sample1.txt")


#for the use of counting number of gene
sumenvarp<-aggregate(table$envarp, by=list(Category=table$gene), FUN=sum)
sumenvarpfc<-aggregate(table$envarpfc, by=list(Category=table$gene), FUN=sum)[,2]
tablesum<-data.frame(cbind(sumenvarp,sumenvarpfc))
colnames(tablesum)<-c("gene","sumenvarp","sumenvarpfc")



hiernormalinvg<-"
data{ #get the data set 
int<lower=0> N;   # number of exon level
int<lower=0> J;   # number of gene level
int <lower=1,upper=J> gene[N];  #index of gene
int <lower=1,upper=N> exon[N];  #index of exon 
vector[N] xij;   #x at exon level
vector[N] yij; #y at exon level
}
parameters{ #specify the parameter we want to know 
real beta;  #common slope for the exon level
real mu;      #common intercept for the exon level
vector[N] aij; #random intercept for the exon level
real <lower=0> sigma_aj2[J];  #variance of intercept at exon level 
vector[J] aj; #random intercept for the gene level 
real <lower=0> sigma_a;  #variance of intercept at gene level
real <lower=0> sigma; #variance of yij
real <lower=0> eps; #hyperparameter for sigma_aj
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
   eps ~ uniform(0,100);
   sigma_aj2 ~inv_gamma(eps,eps);
   aj ~ normal(0, sigma_a);
   for (i in 1:N)
       aij[i]~ normal(0,sqrt(sigma_aj2[gene[i]]));
   yij ~ normal(lambdaij,sigma);
   }

"
library("rstan")
#data 
J<-dim(tablesum)[1] #gene number 
N<-dim(table)[1]  #exon number 
xij=c(table$envarp)
yij=table$envarpfc
gene<-as.numeric(table$gene) #list
genelevel<-length(unique(gene)) #number
indexg<-match(gene, unique(gene))  #list
exon<-c(1:length(table$envarpfc))
M1_table<-list(N=N, J=J, xij=xij,
yij=yij,gene=indexg, exon=exon)
control=list(adapt_delta=0.99,max_treedepth=12)



#3 kinds of iter and warmup 
iter=5000
warmup=1000

fitinv0<-stan(model_code=hiernormalinvg, data=M1_table,iter=iter,warmup=warmup,chains=4) 
print(fitinv0,par='aij')
print(fitinv0,par='aj')
print(fitinv0,par='sigma_aj2')
print(fitinv0,par='sigma_a')
print(fitinv0,par='sigma')
print(fitinv0,par='mu')
print(fitinv0,par='beta')
print(fitinv0,par='eps')

iter=10000
warmup=4000

fitinv1<-stan(model_code=hiernormalinvg, data=M1_table,iter=iter,warmup=warmup,chains=4) 
print(fitinv1,par='aij')
print(fitinv1,par='aj')
print(fitinv1,par='sigma_aj2')
print(fitinv1,par='sigma_a')
print(fitinv1,par='sigma')
print(fitinv1,par='mu')
print(fitinv1,par='beta')
print(fitinv1,par='eps')

iter=20000
warmup=10000

fitinv2<-stan(model_code=hiernormalinvg, data=M1_table,iter=iter,warmup=warmup,chains=4) 
print(fitinv2,par='aij')
print(fitinv2,par='aj')
print(fitinv2,par='sigma_aj2')
print(fitinv2,par='sigma_a')
print(fitinv2,par='sigma')
print(fitinv2,par='mu')
print(fitinv2,par='beta')
print(fitinv2,par='eps')



##############
## end
##############