rm( list = ls()) #clear env
library(escapeMR)
library(dplyr)
library(mra)
#chList <- askForData()
#saveRDS(chList,'chList.rds')
chList<-readRDS('chList.rds')
ch <- ch.for.chat <- chList$ch
covars.exist <- chList$covars.exist
sex=chList$sex
length=chList$length
intervals <- chList$intervals

chat=-1
ns = dim(ch)[2]
nan = dim(ch)[1]
Chat_user = chat

cap.model = as.formula("~ivar(sex,ns)")
surv.model =as.formula("~ivar(length,ns)")

fit_model = Nesc = Bstarall = S = survar1 = survar2 = capvar1 = capvar2 = NULL

if(is.null(intervals)){
  intervals <- rep(1,ncol(ch)-1)
}

if(covars.exist == "YES"){
  fit_model = F.cjs.estim(cap.model, surv.model, histories=ch, c.hat=Chat_user, intervals=intervals)
} 

# nj is the total number of carcasses captured (and checked for marks) at sampling occasion j.
nj = apply(ch,2,function(x){sum(x > 0)})

# chopped are the number amnong the nj that are chopped (no longer part of the marked or unmarked population)
chopped = apply(ch,2,function(x) {sum(x > 1)})

# rj is the total number of carcasses at occasion j that are released with marks. 
rj = nj - chopped

# S is the number of sampling occasions.
Samp = dim(ch)[2]

Bstarall = NULL
for(k in 2:(Samp-2) ) {
  Bj = fit_model$n.hat[k+1] - mean(fit_model$s.hat[,k])*(fit_model$n.hat[k] - (nj[k] - rj[k]) )
  Bstarj = Bj*log(mean(fit_model$s.hat[,k]))/(mean(fit_model$s.hat[,k]) - 1)
  Bstarall = c(Bstarall,Bstarj)
}
Nesc = ceiling(fit_model$n.hat[2]*log(mean(fit_model$s.hat[,1]))/(mean(fit_model$s.hat[,1]) - 1 ) + sum(Bstarall,na.rm=T))



S = dim(ch)[2]

lenInModel <- 0
length.orig <- length
sex.orig <- sex
S = dim(ch)[2]


###
#p.hat and s.hat variance exploring
p.hat<-seq(from=0.05,to=0.95,by=0.05)
s.hat<-seq(from=0.05,to=0.95,by=0.05)
nj<-fit_model$num.caught
rj<-nj - chopped

test_results<-c()
for (i in 1:length(p.hat)){
  for (j in 1:length(s.hat)){
    p<-p.hat[i]
    names(p)='p'
    s<-s.hat[j]
    names(s)='s'
    nhat<-c()
    for(n in 1:S){
      nhat[n]<-nj[n]*(1/p)
    }
    names(nhat)<-c(paste('nhat',1:S,sep=''))
    bhat<-c()
    bstar<-c()
    for(n in 2:S-2){
      bhat[n]<-nhat[n+1]-s*(nhat[n]-(nj[n]-rj[n]))
      bstar[n]<-bhat[n]*log(s)/(s-1)
    }
    names(bhat)<-c(paste('bhat',1:(S-2),sep=''))
    names(bstar)<-c(paste('bstar',1:(S-2),sep=''))
    bstarsum<-sum(bstar)
    names(bstarsum)<-'bstarsum'
    Nesc<-nhat[2]*(log(s)/(s-1))+bstarsum
    names(Nesc)<-'Nesc'
    t<-c(p,s,nhat,bhat,bstar,bstarsum,Nesc)
    test_results<-test_results%>%rbind(t)
  }
}
