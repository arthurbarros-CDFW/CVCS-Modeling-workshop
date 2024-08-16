#speeding up cjs_loglik
#TESTING RCPP FUNCTIONS
library(optimParallel)

source('scripts/function_holder.R')

ch<-read.csv('data/demo_CH.csv')
covars<-read.csv('data/demo_Covariates.csv')
chops<-read.csv('data/demo_Chops.csv')

prepped_data<-CJS_data_prep(ch,chops,covars)

ch_prepped<-prepped_data$ch
cap_X_prepped<-prepped_data$lengths_matrix
surv_X_prepped<-prepped_data$sex_matrix

#testing CJS_loglik

#Get locations outside function

CJS_loglik_testing<-function(beta,ch,cap_X,surv_X){
  #Part 1: Initialize variables
  xlnlik <- 0
  nan=nrow(ch)
  ns=ncol(ch)
  
  first<-location(nan,ns,ch)$first
  last<-location(nan,ns,ch)$last
  
  #Part 3: Calculating total log-likelihood
  for(i in 1:nan){
    
    #set initial values
    sum1=0
    sum2=0
    vp_ij <- matrix(NA, nrow = nan, ncol = ns)
    vs_ij <- matrix(NA, nrow = nan, ncol = ns - 1)
    
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
  }
  return(xlnlik)
}

optim_beta<-c(0.999105767, -0.004392236,  1.431106484, -0.806127998)
CJS_loglik_testing(beta=optim_beta,ch_prepped,cap_X_prepped,surv_X_prepped)

CJS_loglik_wrapper_testing <- function(beta,ch, cap_X, surv_X) {
  -CJS_loglik_testing(beta,ch, cap_X, surv_X)  # Return the negative log-likelihood
}

initial_beta <- numeric(4)

start<-Sys.time()
optim_results <- optim(
  par = initial_beta,
  fn = CJS_loglik_wrapper_testing,
  cap_X = cap_X_prepped,
  surv_X = surv_X_prepped,
  ch = ch_prepped,
  method = 'BFGS',
  control = list(maxit = 1000)
)
end<-Sys.time()
end-start

cl <- makeCluster(detectCores()-1)
setDefaultCluster(cl=cl)# set 'cl' as default cluster

parallel::clusterExport(cl, c("CJS_loglik_testing","location","pro_capsur"))

start<-Sys.time()
optim_results_new <- optimParallel(
  par = initial_beta,
  fn = CJS_loglik_wrapper_testing,
  cap_X = cap_X_prepped,
  surv_X = surv_X_prepped,
  ch = ch_prepped,
  method = 'BFGS',
  lower=c(-Inf, .0001)
)
end<-Sys.time()
end-start #~10 seconds


start<-Sys.time()
nlm(f=CJS_loglik_wrapper_testing,
    p=initial_beta,
    cap_X = cap_X_prepped,
    surv_X = surv_X_prepped,
    ch = ch_prepped
    )
end<-Sys.time()
end-start #~10 seconds