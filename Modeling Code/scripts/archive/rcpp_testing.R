#TESTING RCPP FUNCTIONS
library(Rcpp)
library(RcppArmadillo)
library(roptim)

source('scripts/function_holder.R')

ch<-read.csv('data/demo_CH.csv')
covars<-read.csv('data/demo_Covariates.csv')
chops<-read.csv('data/demo_Chops.csv')

prepped_data<-CJS_data_prep(ch,chops,covars)

ch_prepped<-prepped_data$ch
cap_X_prepped<-prepped_data$lengths_matrix
surv_X_prepped<-prepped_data$sex_matrix

initial_beta<-numeric(4)
optim_beta<-c(0.999105767, -0.004392236,  1.431106484, -0.806127998)


CJS_loglik(optim_beta,cap_X_prepped,surv_X_prepped,ch_prepped) #89.87605

#sourceCpp('Rcpp sources/CJS_functions.cpp')
sourceCpp('Rcpp sources/CJS_functions.cpp')

cpp_location(nan,ns,as.matrix(ch_prepped))
location(nan,ns,as.matrix(ch_prepped))

cpp_pro_capsur(i=3,j=1,beta=optim_beta,as.matrix(ch_prepped),as.matrix(cap_X_prepped),surv_X=as.matrix(surv_X_prepped))
pro_capsur(i=3,j=1,beta=optim_beta,as.matrix(ch_prepped),as.matrix(cap_X_prepped),surv_X=as.matrix(surv_X_prepped))

sourceCpp('Rcpp sources/CJS_functions.cpp')
{gc()
cpp_CJSloglik(optim_beta,as.matrix(ch_prepped),as.matrix(cap_X_prepped),surv_X=as.matrix(surv_X_prepped))
} #this seems to solve a memory leak?


sourceCpp('Rcpp sources/CJS_functions_testing.cpp')
{gc()
  cpp_CJSloglik_testing(optim_beta,as.matrix(ch_prepped),as.matrix(cap_X_prepped),surv_X=as.matrix(surv_X_prepped))
} #trying to deal

###################################################################################
#OPTIM TRIALS
###################################################################################
{gc()
start<-Sys.time()
results1<-optim(par=initial_beta,
              fn=CJS_loglik_wrapper,
              method="BFGS",
              ch=as.matrix(ch_prepped),
              cap_X=as.matrix(cap_X_prepped),
              surv_X=as.matrix(surv_X_prepped),
              control=list(maxit=1000))
end<-Sys.time()

run_time1<-end-start #4 seconds
}

{gc() #often failing?
  start<-Sys.time()
  results2<-optim(par=initial_beta,
                  fn=cpp_CJSwrapper,
                  method="BFGS",
                  ch=as.matrix(ch_prepped),
                  cap_X=as.matrix(cap_X_prepped),
                  surv_X=as.matrix(surv_X_prepped),
                  control=list(maxit=1000))
  end<-Sys.time()
  
  run_time2<-end-start #~8 seconds
}

{gc()
  cpp_optim_testing(beta=initial_beta,
          ch=as.matrix(ch_prepped),
          cap_X=as.matrix(cap_X_prepped),
          surv_X=as.matrix(surv_X_prepped))
}


###################################################################################
#NLMINB TRIALS
###################################################################################

{gc()
start<-Sys.time()
results1<-nlminb(objective=CJS_loglik_wrapper,
                 start=initial_beta,
                 hessian=TRUE,
                 ch=as.matrix(ch_prepped),
                 cap_X=as.matrix(cap_X_prepped),
                 surv_X=as.matrix(surv_X_prepped),
                 lower = -Inf, upper = Inf,
                 control=list(eval.max=1000)
)
end<-Sys.time()
  
run_time1<-end-start #~3
}

#try using cpp_CJS_wrapper
{gc()
  start<-Sys.time()
  results2<-nlminb(objective=cpp_CJSwrapper,
                   start=initial_beta,
                   hessian=TRUE,
                   ch=as.matrix(ch_prepped),
                   cap_X=as.matrix(cap_X_prepped),
                   surv_X=as.matrix(surv_X_prepped),
                   lower = -Inf, upper = Inf,
                   control=list(eval.max=1000)
  )
  end<-Sys.time()
  
  run_time2<-end-start #~3
}

###################################################################################
#BOOTSTRAPPING TRIALS
###################################################################################

catch_matrices<-c()
results<-data.frame()
initial_beta=numeric(4)
iterations<-100

iter_start<-Sys.time()
for(r in 1:iterations){
  iter_start<-Sys.time()
  index=1:dim(ch_prepped)[1]
  samp<- sample(index, replace = T)
  
  ch_iteration = ch_prepped[samp,]
  cap_X_iteration=cap_X_prepped[samp,]
  surv_X_iteration=surv_X_prepped[samp,]
  
  {gc()
    optim_results<-cpp_optim_testing(beta=initial_beta,
              ch=as.matrix(ch_prepped),
              cap_X=as.matrix(cap_X_prepped),
              surv_X=as.matrix(surv_X_prepped))
  }
  
  iter_beta<-optim_results$par
  iter_lik<-optim_results$value
  
  nan=nrow(ch_prepped)
  ns=ncol(ch_prepped)
  p_hat<-s_hat<-matrix( 0, nan, ns )
  
  for(i in 1:nan){
    for(j in 1:ns){
      p_hat[i,j]<-pro_capsur(i,j,ch_iteration,iter_beta,cap_X_iteration,surv_X_iteration)$p.hat
      s_hat[i,j]<-pro_capsur(i,j,ch_iteration,iter_beta,cap_X_iteration,surv_X_iteration)$s.hat
    }
  }
  
  p_hat[,1]<-NA
  
  
  est_escapement<-ceiling(total_escapement(ch_iteration,s_hat,p_hat))
  
  iter_end<-Sys.time()
  iter_time<-iter_end-iter_start
  
  catch_matrices[[r]]<-ch_iteration
  d<-data.frame("iteration"=r,"log-likelihood"=iter_lik,"escapement"=est_escapement,"time"=iter_time)
  results<-results%>%rbind(d)
  print(paste('iteration ',r))
}

hist(results$escapement)
