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
   sigma_aj2 ~inv_gamma(0.05,0.05);
   aj ~ normal(0, sigma_a);
   for (i in 1:N)
       aij[i]~ normal(0,sqrt(sigma_aj2[gene[i]]));
   yij ~ normal(lambdaij,sigma); #log normal 
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
#fitinv<-stan(model_code=hiernormalinvg, data=M1_table,iter=100,chains=1) #10,000 samplings 
print(fitinv)

#plot the density of parameters
#sigma_aj,sigma_a,sigma^2,eps
answer<-extract(fitinv,permuted=TRUE)
print(answer)
print(answer$aij)

plotdes<-function(J,N){
	pdf(file = "inverse gamma (0.5,0.5) prior variance density plot.pdf")
	for (i in 1:J){
		plot(density(answer$sigma_aj2[,i]),main=c("density plot of exon-level variance",i))
		}
	for (i in 1:J){
        plot(density(answer$aj[,i]),main=c("density plot of gene-level intercept",i))
		}
	for (j in 1:N){
		plot(density(answer$aij[,j]),main=c("density plot of exon-level intercept",j))
		}
    plot(density(answer$beta),main="density plot of beta")
	plot(density(answer$mu),main="density plot of mu")	
	plot(density(answer$sigma_a),main="density plot of sigma_a")
	plot(density(answer$sigma),main="density plot of sigma")
	dev.off()
    
}
plotdes(J,N)

