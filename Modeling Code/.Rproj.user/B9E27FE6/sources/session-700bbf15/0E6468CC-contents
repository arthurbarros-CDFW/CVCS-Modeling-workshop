#Methods comparisons
ch<-read.csv('data/additional_data/SRwinter_2021_CH.csv')
covars<-read.csv('data/additional_data/SRwinter_2021_Covariates.csv')
chops<-read.csv('data/additional_data/SRwinter_2021_Chops.csv')

prepped_data<-CJS_data_prep(ch,chops,covars)

ch_prepped<-prepped_data$ch
cap_X_prepped<-prepped_data$lengths_matrix
surv_X_prepped<-prepped_data$sex_matrix


#############################################
#R_optim
#############################################

starttime<-Sys.time()
{gc()
  optim_results<-cpp_optim(beta=initial_beta,
                           ch=as.matrix(ch_prepped),
                           cap_X=as.matrix(cap_X_prepped),
                           surv_X=as.matrix(surv_X_prepped))
}
endtime<-Sys.time()
optim_speed<-endtime-starttime 

est_escapement<-total_escapement(ch_prepped,optim_results$par,
                                 cap_X_prepped,surv_X_prepped)

starttime<-Sys.time()
bootstrap_results<-CJS_bootstrap(iterations=500,
                                 ch_prepped,
                                 cap_X_prepped,
                                 surv_X_prepped)
endtime<-Sys.time()
bootstrap_speed<-endtime-starttime  #7.06 hours for 500 iterations
saveRDS(bootstrap_results,'outputs/SRwinter_2021_bootstrap_results.rds')
bootstrap_results<-readRDS('outputs/SRwinter_2021_bootstrap_results.rds')

#############################################
#escapeMR
#############################################

mra_results<-data.frame()

#Run the for loop for each iteration
for(r in 1:500){
  iter_start<-Sys.time()
  
  #Part 2.1: index and sample the capture histories
  index=1:dim(ch_prepped)[1]
  samp<- sample(index, replace = T)
  
  #Part 2.2: create new capture histories and
  #cap_X and surv_X matrices based on sampled indices
  ch_iteration = ch_prepped[samp,]
  cap_X_iteration=cap_X_prepped[samp,]
  surv_X_iteration=surv_X_prepped[samp,]
  
  #Part 2.3: use cpp_optim to find the optimizal beta
  #parameters for the given iteration data
  source("scripts/escapeMR_testing.R")
  iter_beta<-ans$parameters
  iter_lik<-ans$loglik
  
  #Part 2.4: estimate total escapement for iteration
  est_esc_iter<-total_escapement(ch_iteration,
                                 iter_beta,
                                 cap_X_iteration,
                                 surv_X_iteration)
  
  iter_end<-Sys.time()
  iter_time<-iter_end-iter_start
  
  #Part 2.5: store relevant iteration information in results
  d<-data.frame("iteration"=r,
                "log_likelihood"=iter_lik,
                "escapement"=est_esc_iter,
                "time"=iter_time)
  mra_results<-mra_results%>%rbind(d)
  print(r)
}
saveRDS(mra_results,'outputs/mra_results_SRwinter_2021.rds')
mra_results<-readRDS('outputs/mra_results_SRwinter_2021.rds')
na_results<-mra_results[!complete.cases(mra_results),]