#capture history simulation
#this script will generate simulated capture histories for mark recapture CJS analysis
ch_sim<-function(N1,ns,coef){
  library(tidyverse)
  library(dplyr)
  cap_beta<-coef[1:2]
  surv_beta<-coef[3:4]
  
  #create covariate matrices
  lengths_matrix<-matrix(round(sample(rnorm(100,mean=750,sd=25))),N1,ns)
  sex_matrix<-matrix(sample(c('M','F'),10, prob=c(.5, .5), replace=TRUE),N1,ns)
  sex_matrix<-ifelse(sex_matrix=='F',1,0)
  #survival prediction
  s_model<-function(surv_beta){
    zs<-exp(surv_beta[1]*1+surv_beta[2]*sex_matrix[i,j])
    s.hat<-zs/(1+zs)
  }
  
  p_model<-function(cap_beta){
    zp<-exp(cap_beta[1]*1+cap_beta[2]*lengths_matrix[i,j])
    p.hat<-zp/(1+zp)
  }
  
  p_matrix<-matrix(NA,N1,ns)
  for(i in 1:N1){
    for(j in 1:ns){
      p_matrix[i,j]<-p_model(cap_beta)
    }
  }
  
  s_matrix<-matrix(NA,N1,ns)
  for(i in 1:N1){
    for(j in 1:ns){
      s_matrix[i,j]<-s_model(surv_beta)
    }
  }
  
  #set arrival date
  arrival<-round(rnorm(N1,mean = ns/2,sd=ns/5))
  arrival<-ifelse(arrival>=ns,ns-1,arrival)
  arrival<-ifelse(arrival<1,1,arrival)
  #hist(arrival)
  
  ch <- matrix(NA, N1, ns)
  alive<-matrix(NA, N1, ns)
  alive[,1]=1
  for(j in 1:ns){
    for(i in 1:N1){
      if(arrival[i]>=(j)){ #if fish hasn't arrived, it's still 'alive'
        ch[i,j]=0
        alive[i,j]=1
      } else {
        if(alive[i,(j-1)]==0){ #if dead last period, decide if carcass chopped if wasn't alive at last obs period
          alive[i,j]=0
          ch[i,j]=rbinom(1,2,p_matrix[i,j]*s_matrix[i,j])
        }else{ #determine if fish 'alive', and determine if it was captured or not
          alive[i,j]<-rbinom(1,1,s_matrix[i,j]^sum(alive[i,(arrival[i]:j)],na.rm = T))
          ch[i,j]<-rbinom(1,1,p_matrix[i,j]*s_matrix[i,j])
        }
      }
    }
  }
  for(i in 1:N1){
    for(j in 2:ns){
      if(ch[i,j-1]==2){#if carcass chopped at period j-1, all further values =0
        ch[i,(j:ns)]=0
      }
    }
  }
  
  ID<-1:nrow(ch)
  ch<-ch%>%cbind(ID)
  lengths_matrix<-lengths_matrix%>%cbind(ID)
  sex_matrix<-sex_matrix%>%cbind(ID)
  
  ch<-ch[rowSums(ch[,1:ns])>0,]
  caughtID<-ch[,'ID',drop=F]
  sex_matrix<-subset(sex_matrix, ID %in% caughtID)
  lengths_matrix<-subset(lengths_matrix, ID %in% caughtID)
  
  ch<-ch[,-(ns+1)]
  sex_matrix<-sex_matrix[,-(ns+1)]
  lengths_matrix<-lengths_matrix[,-(ns+1)]
  
  sim_list<-list('ch'=ch,'sex'=sex_matrix,'lengths'=lengths_matrix,'arrivals'=arrival,'alive'=alive)
  return(sim_list)
}

