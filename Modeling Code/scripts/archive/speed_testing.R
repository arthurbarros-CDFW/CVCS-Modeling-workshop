#Speed tests for optimization algorithms
#extended model comparisons
rm( list = ls()) #clear env
library(dplyr)
library(stats4)
library(escapeMR)
library(ggplot2)
library(mra)
library(Rcpp)
library(roptim)

source('scripts/ch_simulator.R')
source('scripts/CJS_functions.R')
source('scripts/schnabel_est.R')
sourceCpp('scripts/Rcpp sources/CJS_functions.cpp')

############################################################################
#Simulation testing
############################################################################
sim_coef<-c(3.908292515, -0.006356866,  2.972067636, -0.105286686)
N1_range<-seq(from=1000,to=5000,by=1000)
ns_range<-seq(from=15,to=45,by=10)
library(memoise)
#apparently this memoise package speeds things up, not sure how it works but it does, note that it only speeds up after the initial time the function is called
data<-NULL
schnabel_results<-NULL

for(q in 1:length(N1_range)){
  for(w in 1:length(ns_range)){
    N1=N1_range[q]
    ns=ns_range[w]
    iterations<-1
    results<-data.frame()
    iter_start<-Sys.time()
    catch_matrices<-c()
    sim<-ch_sim(N1,ns,sim_coef)
    catch_hist<-data.frame(sim$ch)
    sex_matrix<-sim$sex
    lengths_matrix<-sim$lengths
    DiscTag<-rep(1:nrow(catch_hist))
    covariates<-data.frame(DiscTag,sex=sex_matrix[,1],lengths=lengths_matrix[,1])
    write_catch<-data.frame(DiscTag,catch_hist)
    mf<-memoise(CJS_obj)
    write.csv(write_catch,'outputs/sim_data/sim_ch.csv',row.names = F)
    write.csv(covariates,'outputs/sim_data/sim_covariates.csv',row.names = F)
    
    for(r in 1:iterations){
      iter<-r
      index = 1:dim(catch_hist)[1]
      samp<- sample(index, replace = T)
      ch = catch_hist[samp,]
      
      #Set Global Variables, can't figure out why I can't pass these in the function...
      nan=nrow(ch)
      ns=ncol(ch)
      p=np=4
      nx=ny=p/2
      beta= numeric(p)
      
      cap_X=as.matrix(lengths_matrix[samp,])
      surv_X=as.matrix(sex_matrix[samp,])
      
      #R-script CJS
      print(c(paste('N1=',N1),paste('ns=',ns),paste('Total Marked=',nrow(ch))))
      lengths_vector<-cap_X[,1]
      sex_vector<-surv_X[,1]
      #but let's also test all the different maximization methods available in optim
      om<-c("R CJS")
      iter_df<-data.frame()
      starttime<-Sys.time()
      CJS_coef<-tryCatch({optim(beta, CJS_obj, method = 'BFGS',hessian=T)$par
      },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
      )
      if(is.null(CJS_coef)){
        CJStime<-CJS_lik<-CJS_N<-'error'
      }else{
        CJS_N<-CJS_escapement(ch,CJS_coef,np)
        CJS_lik<-CJS_loglik(CJS_coef,p)
        endtime<-Sys.time()
        CJStime<-endtime-starttime
      }
      d<-data.frame(N1,iter,om,CJS_N,CJS_lik,CJStime)
      iter_df<-iter_df%>%rbind(d)
      ch_matrix = data.matrix(ch)
      
      #RCPP CJS
      starttime<-Sys.time()
      om<-c("RCPP CJS")
      CJS_coef<-tryCatch({cpp_optim(
        beta=beta,
        ns=ns,
        nan=nan,
        ch=ch_matrix,
        p=p,
        cap_X=cap_X,
        surv_X=surv_X,
        nx=nx,
        np=np
      )$par},error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
      )
      endtime<-Sys.time()
      CJStime<-endtime-starttime
      
      if(is.null(CJS_coef)){
        CJStime<-CJS_lik<-CJS_N<-'error'
      }else{
        CJS_N<-CJS_escapement(ch,CJS_coef,np)
        CJS_lik<-CJS_loglik(CJS_coef,p)
        CJStime<-endtime-starttime
      }
      d<-data.frame(N1,iter,om,CJS_N,CJS_lik,CJStime)
      iter_df<-iter_df%>%rbind(d)
      
      #F.cjs_estim
      starttime<-Sys.time()
      source('scripts/package_testing/F.cjs.estim_testing.R')
      CJS_coef<-ans$parameters
      CJS_N<-CJS_escapement(ch,CJS_coef,np)
      CJS_lik<-ans$loglik
      endtime<-Sys.time()
      CJStime<-endtime-starttime
      d<-data.frame(N1,iter,om='MRA',CJS_N,CJS_lik,CJStime)
      iter_df<-iter_df%>%rbind(d)
      
      #get schnabel estimate
      N_schn<-schnabel(ch,ns,nan)[1]
      CI_schn<-schnabel(ch,ns,nan)[2]
      schn_data<-data.frame(N_schn,CI_schn)
      
      #iter_df<-iter_df%>%cbind(schn_data)
      results<-results%>%rbind(iter_df)
    }
    results$observations<-ns
    results$total_marked=nrow(ch)
    data<-data%>%rbind(results)
    saveRDS(data,'outputs/test_sim_results.rds')
    drop_cache(mf)
  }
}

