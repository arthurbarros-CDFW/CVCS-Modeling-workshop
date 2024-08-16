library(knitr)
library(tinytex)
library(tidyverse)

CJS_loglik<-function(beta,cap_X,surv_X,ch){
  row_sums<-data.frame()
  #Part 1: Initialize variables
  xlnlik <- 0
  nan=nrow(ch)
  ns=ncol(ch)
  
  #Part 2: get locations
  first<-location(nan,ns,ch)$first
  last<-location(nan,ns,ch)$last
  
  #Part 3: Calculating total log-likelihood
  for(i in 1:nan){
    
    #set initial values
    sum1=0
    sum2=0
    vp_ij<-NA
    vs_ij<-NA
    
    if(first[i]==0){
      init_cap=ns+1
      init_surv=ns+1
    } else if(first[i]>0){
      init_cap=first[i]
      init_surv=first[i]-1
    }
    
    #Part 3.1: Compute all probabilities of capturing carcass i,
    #from first occasion to end
    if (init_cap<=ns){
      for (j in init_cap:ns) {
        vp_ij[j] <- pro_capsur(i,j,ch,beta,cap_X,surv_X)$p.hat
      }
    }
    #Part 3.2:Compute all probabilities of carcass i surviving from j to j+1,
    #from first occasion to end
    if(init_surv<ns){
      for(j in init_surv:(ns-1)){
        vs_ij[j] <- pro_capsur(i,j,ch,beta,cap_X,surv_X)$s.hat
      }
    }
    
    #Part 3.3: compute log-likelihood contribution for animal i
    if((first[i]>0 && first[i]<=last[i])){
      for(j in first[i]:last[i]){#for each of the observations for the animal,
        #from first contact to last
        hij <- ifelse(ch[i, j] >= 1, 1, 0)
        #for each observation, calculate a running sum of
        #ln(L)=ln(p)+(ln(1-p))+ln(survival)
        sum1=sum1+hij*log(vp_ij[j])+
          (1-hij)*log(1-vp_ij[j])+
          log(vs_ij[j-1])
      }
    }
    
    #Part 3.4: Find second part of likelihood, Chi
    #Chi = probability of animal i not being seen again after last i
    #If animal died on capture before release, prob of not seeing again is 1
    if(ch[i,last[i]]>=2){
      sum2=0
    }else if(last[i]>0 && last[i]<ns){
      sum2=1-vs_ij[last[i]] #chance it died at last[i]
      for(ii in (last[i]+1):ns){#for each obs after last[i]
        prod<-1
        for(jj in last[i]:(ii-1)){
          prod<-prod*vs_ij[jj]*(1-vp_ij[jj+1])#product of all chances carcass
          #survived but wasn't seen for each observation point between
          #last and ns-1
        }
        if(ii<ns){
          prod<-prod*(1-vs_ij[ii])
        }
        sum2<-sum2+prod
      }
      sum2<-log(sum2)
    }
    
    #Part 3.5: estimate total likelihood of carcass capture history
    xlnlik<-xlnlik+sum1+sum2
    print(c(paste("i: ",i),paste("sum1: ",sum1),paste("sum2: ",sum2),paste("xlnlik: ",xlnlik)))
  }
  return(xlnlik)
}

pro_capsur<-function(i,j,ch, beta,cap_X,surv_X){
  nan=nrow(ch)
  ns=ncol(ch)
  
  p=length(beta)
  
  #purpose: evaluate probability of capture and survival for each animal i
  cap_beta<-beta[1:(p/2)]
  surv_beta<-beta[((p/2)+1):p]
  
  zp<-exp(cap_beta[1]*1+cap_beta[2]*cap_X[i,j])
  zs<-exp(surv_beta[1]*1+surv_beta[2]*surv_X[i,j])
  
  p.hat<-zp/(1+zp)
  s.hat<-zs/(1+zs)
  
  est_list<-list('p.hat'=p.hat,'s.hat'=s.hat)
  return(est_list)
}

location<-function(nan,ns,ch){
  #purpose: compute first and last capture for each animal
  first <- rep(0, nan)
  last <- rep(0, nan)
  for(i in 1:nan){
    findch=TRUE
    for(j in 1:ns){
      if(ch[i,j]>=1){
        if(findch==T && (j<ns)){
          first[i]=(j+1)
          findch=FALSE
        }
        last[i]=j
      }
    }
  }
  est_list<-list('first'=first,'last'=last)
  return(est_list)
}

CJS_data_prep<-function(ch,chops,covars){
  set.seed(4211)
  #Part 1: prepare ch and covars data
  ch=ch[-1] #remove disctag vector from capture histories
  covars$sex<-ifelse(covars$sex=='F',1,0) #change sex to numeric value 
  
  #Part 2: prep chops data
  chops<-chops[-1]
  clean_chops<-matrix(ncol=ncol(chops))
  for(i in 1:ncol(chops)){
    d<-chops[i]
    r<-rep(0,ncol(chops))
    r[i]=2
    r<-matrix(rep((r),d),ncol=ncol(chops),byrow=T)
    clean_chops<-clean_chops%>%rbind(r)
  }
  clean_chops<-clean_chops[-1,]
  colnames(clean_chops)<-colnames(ch)
  #add in chops
  ch<-ch%>%rbind(clean_chops)
  
  #Part 3: generate sex data for chopped data
  index = 1:dim(ch)[1]
  samp<- sample(index, replace = T)
  
  covariates_new<-as.matrix(covars[samp,])
  
  sex_vector<-covariates_new[,'sex']
  no.miss = sex_vector[!is.na(sex_vector)]
  prop.female<-mean(no.miss)
  prop.male<-1-prop.female
  sex_vector<-ifelse(is.na(sex_vector),sample(c(1,0),1,prob=c(prop.female,prop.male)),sex_vector)
  sex_matrix<-matrix(sex_vector,nrow=nrow(ch),ncol=ncol(ch))
  
  #Part 4: generate length data for chopped data
  lengths_vector<-covariates_new[,'length']
  females<-covariates_new[covariates_new[,2]==1,]
  males<-covariates_new[covariates_new[,2]==0,]
  avg.female.length<-round(mean(as.numeric(females[,3]),na.rm=T),1)
  avg.male.length<-round(mean(as.numeric(males[,3]),na.rm=T),1)
  for(i in 1:length(lengths_vector)){
    if(is.na(lengths_vector[i])){
      lengths_vector[i] = ifelse(sex_matrix[i,1] == 1,  avg.female.length, avg.male.length)
    }
  }
  lengths_matrix<-matrix(as.numeric(lengths_vector),nrow=nrow(ch),ncol=ncol(ch))
  
  return(list('ch'=ch,'lengths_matrix'=lengths_matrix,'sex_matrix'=sex_matrix))
}

CJS_loglik_wrapper <- function(beta, cap_X, surv_X, ch) {
  -CJS_loglik(beta, cap_X, surv_X, ch)  # Return the negative log-likelihood
}

B_star<-function(ch,s_hat,p_hat){
  R=n=list()
  nan=nrow(ch)
  ns=ncol(ch)
  for(j in 1:ns){
    d<-ch #select just the capture matrix data
    R[j]<-as.numeric(length(which(d[j]==1)))
    n[j]<-as.numeric(length(which(d[j]==1))+length(which(d[j]==2)))
  }
  
  N_hat<-Horvitz_Thompson(p_hat,ch)
  #next B1, or total number of births for each period
  B1<-list()
  for(j in 2:(ns-2)){
    B1[j]<-N_hat[[j+1]]-mean(s_hat[,j])*(N_hat[[j]]-(n[[j]]-R[[j]]))
  }
  
  #next Bstar, number of births adjusted for those entering the system between
  #j and j+1, but not surviving to j+1
  Bstar<-NULL
  for(j in 2:(ns-2)){
    Bstar[j]<-as.numeric(B1[[j]]*(log(mean(s_hat[,j]))/(mean(s_hat[,j])-1)))
  }
  Bstar<-Bstar[-(1)]
  return(Bstar)
}

total_escapement<-function(ch,s_hat,p_hat){
  
  N_hat=Horvitz_Thompson(p_hat,ch)
  
  Bstar=B_star(ch,s_hat,p_hat)
  
  escapement<-N_hat[[2]]*(log(mean(s_hat[,1]))/(mean(s_hat[,1])-1)) + 
    sum(Bstar,na.rm=T)
  
  return(escapement)
}

Horvitz_Thompson<-function(p_hat,ch){
  nan=nrow(ch)
  ns=ncol(ch)
  N_hat<-list()
  n_mat<-matrix(NA,nan,ns)
  for(j in 1:ns){
    for(i in 1:nan){
      n_mat[[i,j]]<-if(ch[[i,j]]>=1){
        1/p_hat[[i,j]]
      } else {0}
    }
    N_hat[[j]]=sum(n_mat[,j])
  }
  return(N_hat)
}
