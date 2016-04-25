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
fitinv<-stan(model_code=hiernormalinvg, data=M1_table,iter=40000,warmup=35000,chains=4) #10,000 samplings 
#fitinv<-stan(model_code=hiernormalinvg, data=M1_table,iter=2000,chains=4) #10,000 samplings 
print(fitinv)


answer<-extract(fitinv,permuted=TRUE)#extract samples

#get the mode of the sample distribution
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#creat the mean, median, mode table
mode=rep(0,N)
mean=rep(0,N)
median=rep(0,N)
for ( i in 1: N){
      mode[i]=Mode(answer$aij[,i]+answer$aj[,indexg[i]])
	  mean[i]=mean(answer$aij[,i]+answer$aj[,indexg[i]])
	  median[i]=median(answer$aij[,i]+answer$aj[,indexg[i]])
}
summary=cbind(mode,mean,median,indexg)
print(summary)
write.table (summary, "3msummary.txt",sep="\t")
