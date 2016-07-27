#####################################################################
## Data management
## 0726 Shuaiqi 

setwd("")
table <- read.table("exon_level_process_v2.txt")
colnames(table) <- c("chr", "gene", "dom", "subdom", "exon", "gene.dom", 
             "gene.dom.subdom", 
             "envarp",    # pass
             "envarpf",   # pass functional
             "envarpfr",  # pass functional rare
             "emutr")     # mutation rate 

table <- within(table,envarpfc<-envarpf-envarpfr)#y
table <- within(table,gene<-factor(gene))
table <- within(table,gene.dom<-factor(gene.dom))
table <- within(table,gene.dom.subdom<-factor(gene.dom.subdom))
table <- table[which(table$envarp!=0), ]
table$x = scale(table$envarp)


#function to get the gene factors as a column 
expanddata <- function(data,gene){
	n <- length(unique(gene))

	for (i in 1:n){
		columnname<-unique(gene)[i]
		newcolumn<-rep("0",dim(data)[1])
		position<-which(gene%in%columnname)
		newcolumn[position]<-as.character(data$dom[position])
		data$newcolumn<-newcolumn
		names(data)[names(data)=="newcolumn"] <- gsub("-","",toString(unique(gene)[i]))
				}
data
}

data <- table
gene <- data$gene
expdata <- expanddata(data,gene)
expdata$y <- log(expdata$envarpfc+0.01)
write.table(expdata, file = "expand data.txt")

########################################################################

### lmer model
#################################################
library("lme4")

setwd("")
expdata <- read.table("expand data.txt")
#get the formular
formular<-function(gene){
  n<-length(unique(gene))
  formular<-"y~1+x+(1|gene)"
  for ( i in 1:n){
    group<-gsub("-","",unique(gene)[i])
    formular<-paste(formular,"+",'(1|',group,')')
  }
  formular
}

gene <- expdata$gene
form <- formular(gene)



#change the designe matrix 
#options(expressions=1000000)
lmod <- lFormula(form,data=expdata)

#options(max.print= 10000, width = 100)

#change Zt matrix 
Z <- lmod$reTrms$Zt
name <- dimnames(Z)[[1]] #get the factor names 
position <- which(name%in%"0")
lmod$reTrms$Zt[position,] <- rep(0,dim(Z)[2])#change the sparse matrix for "0" group 


#change Ztlist
Ztlistc <- function(gene,Ztlist,Z){
  n <- length(unique(gene))
  for ( i in 2:n){
    Ztlist[[i]]["0",] <- t(as(rep(0,dim(Z)[2]),Class="sparseMatrix"))
  }
  Ztlist
}

Ztlist <- lmod$reTrms$Ztlist
lmod$reTrms$Ztlist <- Ztlistc(gene,Ztlist,Z)


#get the result 
dfun <- do.call(mkLmerDevfun,lmod)   ## create dev fun from modified lf

#devfunw <- function(theta) {
#   n <- length(lower)  ## from environment
#   th <- numeric(n)
#   diag_el <- which(lower == 0)
#   th[diag_el] <- theta
#   dfun(th)
#}
#environment(devfunw) <- environment(dfun)

opt<-optimizeLmer(dfun,optimizer="Nelder_Mead",control=list(maxfun=5000000))

##  opt <- Nelder_Mead(devfunw,control=list(maxfun=5000000),par=lmod$reTrms$theta,
  ##  lower = rep.int(0,length(lmod$reTrms$theta)), upper = rep.int(Inf, length(lmod$reTrms$theta)) )   

opt$fval  #numeric scalar - the minimum function value achieved
opt$convergence #integer valued scalar, if not 0, an error code


## make the results into a 'merMod' object
fit <- mkMerMod(environment(dfun), opt, lmod$reTrms,
                fr = lmod$fr)

#write result to txt file         
out <- capture.output(fit)

cat(out,file="sdrandomeffect.txt",sep="\n",append=TRUE)

randomgene <- capture.output(ranef(fit))

cat(randomgene, file = "randomintercept.txt",sep="\n",append=TRUE)


print(as.data.frame(VarCorr(fit)))#covaraince matrix between random intercepts 


## END
########################################################################################################

