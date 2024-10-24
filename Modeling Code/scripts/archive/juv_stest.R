library(escapeMR)
ch<-read.csv('data/additional_data/juv data/caphist2020.csv')
covars<-read.csv('data/additional_data/juv data/covars.csv')
chops<-read.csv('data/additional_data/juv data/chops.csv')
prepped_data<-CJS_data_prep(ch,chops,covars)
ch_prepped<-prepped_data$ch
cap_X_prepped<-prepped_data$lengths_matrix
surv_X_prepped<-prepped_data$sex_matrix

initial_beta=numeric(4)
starttime<-Sys.time()
optim_results<-optim(par=initial_beta,
                     fn=CJS_loglik_wrapper,
                     method="BFGS",
                     ch=as.matrix(ch_prepped),
                     cap_X=as.matrix(cap_X_prepped),
                     surv_X=as.matrix((surv_X_prepped)))
endtime<-Sys.time()
optim_speed_R<-endtime-starttime


optim_beta<-optim_results$par


est_list<-fill_prob_matrices(ch_prepped,optim_beta,cap_X_prepped,surv_X_prepped)

est_escapement<-total_escapement(ch_prepped,optim_beta,
                                 cap_X_prepped,surv_X_prepped)
