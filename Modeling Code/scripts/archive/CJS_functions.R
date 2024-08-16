############################################################################
#establish functions that we will use
############################################################################

#note: this depends on how I structure my covariate equations, and will need to be changed to allow variation on regression equation
pro_capsur<-function(beta,i,j,p){
  #purpose: evaluate probability of capture and survival for each animal i
  #global vars: lengths_matrix,sex_matrix,nx,nan,ns
  cap_beta<-beta[1:nx]
  surv_beta<-beta[(nx+1):np]
  zp<-exp(cap_beta[1]*1+cap_beta[2]*cap_X[i,j])
  p.hat<-zp/(1+zp)
  zs<-exp(surv_beta[1]*1+surv_beta[2]*surv_X[i,j])
  s.hat<-zs/(1+zs)
  est_list<-list('p.hat'=p.hat,'s.hat'=s.hat)
  return(est_list)
}

location<-function(nan,ns,ch){
  #purpose: compute first and last capture for each animal
  first <- rep(0, nan)
  last <- rep(0, nan)
  for(i in 1:nan){
    findic=TRUE
    for(j in 1:ns){
      if(ch[i,j]>=1){
        if(findic==T && (j<ns)){
          first[i]=(j+1)
          findic=FALSE
        }
        last[i]=j
      }
    }
  }
  est_list<-list('first'=first,'last'=last)
  return(est_list)
}


CJS_loglik<-function(beta,p){
  #global vars: cap_matrix,surv_matrix,nx,nan,ns, ch
  #get locations
  first<-location(nan,ns,ch)$first
  last<-location(nan,ns,ch)$last
  #initialize init1 and init2 to kill uninitialized warning at compile
  #Assuming first(i) >= 0, these values are always overwritten in do loop below.
  # Initialize variables
  xlnlik <- 0
  init1 = 1
  init2 = 1
  for(i in 1:nan){
    #First part of the log-likelihood function, between first and last 1
    sum1=0
    sum2=0
    vp_ij<-NA
    vs_ij<-NA
    
    if(first[i]==0){
      init1=ns+1
      init2=ns+1
    } else if(first[i]>0){
      init1=first[i]
      init2=first[i]-1
    }
    
    #Compute all probabilities of capturing bird i, from first occasion to end
    if (init1<=ns){
      for (j in init1:ns) {
        vp_ij[j] <- pro_capsur(beta,i, j,p)$p.hat
      }
    }
    #compute probability of bird i surviving to time j
    if(init2<ns){
      for(j in init2:(ns-1)){
        vs_ij[j] <- pro_capsur(beta,i, j,p)$s.hat
      }
    }
    
    #compute log-likelihood contribution for animal i
    if((first[i]>0 && first[i]<=last[i])){
      for(j in first[i]:last[i]){#for each of the observations for the animal, from first contact to last
       hij <- ifelse(ch[i, j] >= 1, 1, 0)
        #for each observation, calculate a running sum of ln(L)=ln(p)+(ln(1-p))+ln(survival)
        sum1=sum1+hij*log(vp_ij[j])+
          (1-hij)*log(1-vp_ij[j])+
          log(vs_ij[j-1])
      }
    }
    #Find second part of likelihood, after last 1 = Chi parameters
    #Chi = probability of animal i not being seen again
    #If animal died on capture before release, prob of not seeing again is 1
    if(ch[i,last[i]]>=2){
      sum2=0
    }else if(last[i]>0 && last[i]<ns){
      sum2=1-vs_ij[last[i]] #chance it died at last[i]
      for(ii in (last[i]+1):ns){#for each obs after last[i]
        prod<-1
        for(jj in last[i]:(ii-1)){
          prod<-prod*vs_ij[jj]*(1-vp_ij[jj+1])#product of all chances fish survived but wasn't seen for each observation point between last and ns-1
        }
        if(ii<ns){
          prod<-prod*(1-vs_ij[ii])
        }
        sum2<-sum2+prod
      }
      sum2<-log(sum2)
    }
    xlnlik<-xlnlik+sum1+sum2
  }
  return(xlnlik)
}

#################################################
CJS_loglik_testing<-function(beta,p){
  #global vars: cap_matrix,surv_matrix,nx,nan,ns, ch
  #get locations
  first<-location(nan,ns,ch)$first
  last<-location(nan,ns,ch)$last
  #initialize init1 and init2 to kill uninitialized warning at compile
  #Assuming first(i) >= 0, these values are always overwritten in do loop below.
  # Initialize variables
  xlnlik <- 0
  init1 = 1
  init2 = 1
  for(i in 1:nan){
    #First part of the log-likelihood function, between first and last 1
    sum1=0
    sum2=0
    vp_ij<-NA
    vs_ij<-NA
    
    # First(i) is never 1, because subroutine location sets first to 1 occasion
    # after the first encounter = first estimable capture probability.
    # Note first(i) = 0 for histories with first capture at the last occasion
    if(first[i]==0){
      init1=ns+1
      init2=ns+1
    } else if(first[i]>0){
      init1=first[i]
      init2=first[i]-1
    }
    
    #Compute all probabilities of capturing bird i, from first occasion to end
    if (init1<=ns){
      for (j in init1:ns) {
        vp_ij[j] <- pro_capsur(beta,i, j,p)$p.hat
      }
    }
    #compute probability of bird i surviving to time j
    if(init2<ns){
      for(j in init2:(ns-1)){
        vs_ij[j] <- pro_capsur(beta,i, j,p)$s.hat
      }
    }
    
    #compute log-likelihood contribution for animal i
    if((first[i]>0 && first[i]<=last[i])){
      for(j in first[i]:last[i]){#for each of the observations for the animal, from first contact to last
        hij <- ifelse(ch[i, j] >= 1, 1, 0)
        #for each observation, calculate a running sum of ln(L)=ln(p)+(ln(1-p))+ln(survival)
        sum1=sum1+hij*log(vp_ij[j])+
          (1-hij)*log(1-vp_ij[j])+
          log(vs_ij[j-1])
      }
    }
    #Find second part of likelihood, after last 1 = Chi parameters
    #Chi = probability of animal i not being seen again
    #If animal died on capture before release, prob of not seeing again is 1
    if(ch[i,last[i]]>=2){
      sum2=0
    }else if(last[i]>0 && last[i]<ns){
      sum2=1-vs_ij[last[i]] #chance it died at last[i]
      for(ii in (last[i]+1):ns){#for each obs after last[i]
        prod<-1
        for(jj in last[i]:(ii-1)){
          prod<-prod*vs_ij[jj]*(1-vp_ij[jj+1])#product of all chances fish survived but wasn't seen for each observation point between last and ns-1
        }
        
        if(ii<ns){
          prod<-prod*(1-vs_ij[ii])
        }
        
        sum2<-sum2+prod
      }
      sum2<-log(sum2)
      
    }
    xlnlik<-xlnlik+sum1+sum2
  }
    return(xlnlik)
}
#################################################


CJS_obj<-function(beta,p){
  #Purpose: to calculate the Capture-recapture log-likelihood given data and values of parameters.  This is the routine that will be optimized.
  #Calculate the log-likelihood
  lnlik =-1.0 * CJS_loglik(beta,p)
  
  return(lnlik)
  #Calculate the gradient
  #grad <- CJS_gradient(beta, p,lnlik, grad)
  #return(list(lnlik = lnlik, grad = grad))
}

CJS_escapement<-function(ch,beta,np){

  #Take coefficients and produce p.hat/s.hat
  p_hat<-MRA_p_hat<-s_hat<-matrix( 0, nan, ns )
  
  for(i in 1:nan){
    for(j in 1:ns){
      p_hat[i,j]<-pro_capsur(beta,i,j,np)$p.hat
      s_hat[i,j]<-pro_capsur(beta,i,j,np)$s.hat
    }
  }
  #set first row of capture probabilities to NA
  p_hat[,1]<-NA
  
  #Calculate N_hat using Horvitz-Thompson estimator
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
  
  #Estimate births and then escapement
  #calculate R and n
  #R= #the total number of carcasses during j that are released with marks
  #r= #the number of R captured again at some later time
  R=n=list()
  for(j in 1:ns){
    d<-ch #select just the capture matrix data
    R[j]<-as.numeric(length(which(d[j]==1)))
    n[j]<-as.numeric(length(which(d[j]==1))+length(which(d[j]==2)))
  }
  
  #next B1, or total number of births for each period
  B1<-list()
  for(j in 2:(ns-2)){
    B1[j]<-N_hat[[j+1]]-mean(s_hat[,j])*(N_hat[[j]]-(n[[j]]-R[[j]]))
  }
  
  #next B2, number of births adjusted for those entering the system between j and j+1, but not surviving to j+1
  B2<-NULL
  for(j in 2:(ns-2)){
    B2[j]<-as.numeric(B1[[j]]*(log(mean(s_hat[,j]))/(mean(s_hat[,j])-1)))
  }
  B2<-B2[-(1)]
  
  #finally, estimate total escapement
  escapement<-N_hat[[2]]*(log(mean(s_hat[,1]))/(mean(s_hat[,1])-1)) + sum(B2,na.rm=T)
  
  return(escapement)
}

#below not used yet

CJS_gradient <- function(beta,p, f,grad) {
  # Purpose: compute the derivative of the likelihood at beta.
  #
  # Input:
  # p = number of coefficients
  # beta = coefficients
  # f = value of loglikelihood at beta
  #
  # Output:
  # grad = pX1 vector of partial derivatives of log likelihood w.r.t. all parameters.
  #        Derivatives computed by either central differences or one-sided differences.
  #        Method is controlled by the constant central_diffs.  Central differences are
  #        more accurate, but you have to compute the likelihood twice for each parameter,
  #        instead of just once.
  
  grad <- numeric(p)
  beta2 <- beta
  
  deltax <- 1.0e-8
  central_diffs<-TRUE
  
  # central_diffs and deltax set externally
  
  for (i in 1:p) {
    deltai <- (deltax/2.0) * (1.0 + abs(beta2[i])) * 1e5
    
    tmp_b <- beta2[i]
    beta2[i] <- beta2[i] + deltai
    f1 <- -1.0 * CJS_loglik(beta2,p)
    
    if (central_diffs) {
      beta2[i] <- tmp_b - deltai
      f2 <- -1.0 * CJS_loglik(beta2,p)
      grad[i] <- (f1 - f2) / (2.0 * deltai)
    } else {
      grad[i] <- (f1 - f) / deltai
    }
    
    beta2[i] <- tmp_b
  }
  
  return(grad)
}
