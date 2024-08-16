#extended model comparisons
rm( list = ls()) #clear env
library(dplyr)
library(stats4)
library(escapeMR)
library(ggplot2)
library(memoise)
library(mra)

source('scripts/package_testing/ch_simulator.R')
source('scripts/package_testing/CJS_functions.R')
source('scripts/package_testing/archive/capture_simulation.R')
source('scripts/package_testing/schnabel_est.R')
############################################################################
#Simulate population and bootstrap test
############################################################################
MRA_coef<-c(3.908292515, -0.006356866,  2.972067636, -0.105286686)
N1=1000
ns=15
sim_coef<-MRA_coef#just using coefficients produced by MRA for prior data
iterations<-1
results<-data.frame()
iter_start<-Sys.time()
catch_matrices<-c()
library(memoise)
mf<-memoise(CJS_obj) #apparently this memoise package speeds things up, not sure how it works but it does, note that it only speeds up after the initial time the function is called

sim<-ch_sim(N1,ns,sim_coef)
catch_hist<-data.frame(sim$ch)
sex_matrix<-sim$sex
lengths_matrix<-sim$lengths
total_marked<-nrow(catch_hist)


for(r in 1:iterations){
  iter<-r
  print(iter)
  index = 1:dim(catch_hist)[1]
  samp<- sample(index, replace = T)
  ch = catch_hist[samp,]
  
  #Set Global Variables, can't figure out why I can't pass these in the function...
  nan=nrow(ch)
  ns=ncol(ch)
  p=np=4
  nx=ny=p/2
  beta= numeric(p)
  
  sex_orig<-as.matrix(sex_matrix[samp,])
  
  lengths_orig<-as.matrix(lengths_matrix[samp,])
  
  cap_X=lengths_orig
  surv_X=sex_orig
  
  lengths_vector<-cap_X[,1]
  sex_vector<-surv_X[,1]
  #but let's also test all the different maximization methods available in optim
  optim_methods<-c("Nelder-Mead", "BFGS")
  iter_df<-data.frame()
  for(i in 2:length(optim_methods)){
    om<-optim_methods[i]
    starttime<-Sys.time()
    CJS_coef<-tryCatch({optim(beta, CJS_obj, method = om,hessian=T)$par
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
  }
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
  results<-results%>%rbind(iter_df)
  catch_matrices[[r]]<-ch
  drop_cache(mf)
}
iter_end<-Sys.time()
iter_totaltime<-iter_end-iter_start

#confidence intervals
conf.level = 95
alpha = 1 - conf.level/100
lower = alpha/2
upper = 1 - alpha/2
mid=.5

ci<-results
ci<-ci%>%
  group_by(om)%>%
  summarise(lower_ci=ceiling(quantile(CJS_N,probs = c(lower),na.rm=T)),mid_ci=ceiling(quantile(CJS_N,probs = c(mid),na.rm=T)),upper_ci=ceiling(quantile(CJS_N,probs = c(upper),na.rm=T)))

#get schnabel estimate
N_schn<-schnabel(ch,ns,nan)[1]
CI_schn<-schnabel(ch,ns,nan)[2]

#plot results
data<-results
data$CJS_N<-as.numeric(data$CJS_N)
ggplot(data,aes(x=CJS_N))+
  geom_histogram(color = "#000000", fill = "#0099F8")+
  facet_grid(factor(om,levels=c('Nelder-Mead','BFGS','MRA'))~.)+
  geom_segment(data=ci,aes(x=lower_ci,xend=lower_ci,y=0,yend=15),linewidth=1,linetype='dashed')+
  geom_segment(data=ci,aes(x=upper_ci,xend=upper_ci,y=0,yend=15),linewidth=1,linetype='dashed')+
  geom_segment(data=ci,aes(x=mid_ci,xend=mid_ci,y=0,yend=15),linewidth=1,linetype='dashed')+
  geom_vline(xintercept=N1,color='red')+
  geom_text(data=ci,aes(x=mid_ci,y=25,label=paste('Mean:',mid_ci)))+
  geom_text(data=ci,aes(x=lower_ci,y=20,label=paste('Lower:',lower_ci)))+
  geom_text(data=ci,aes(x=upper_ci,y=15,label=paste('Upper:',upper_ci)))+
  theme_classic()
ggsave('figures/sim_comparison_test.png',width=8,height = 5)

############################################################################
#Simulation testing
############################################################################
N1_range<-seq(from=500,to=4000,by=500)
ns_range<-seq(from=5,to=40,by=5)
results<-data.frame()
for(n in 1:length(N1_range)){
  for(o in 1:length(ns_range)){
    d<-data.frame()
    N1=N1_range[n]
    ns=ns_range[o]
    sim<-ch_sim(N1,ns,coef=MRA_coef)
    ch=sim$ch
    lengths=sim$lengths
    sex=sim$sex
    alive<-sim$alive
    dead<-NA
    for(i in 1:N1){
      for(j in 2:ns){
        if(alive[i,j]==0 & alive[i,(j-1)]==1){
          dead[i]<-j
        }
      }
    }
    total_marked=nrow(ch)
    left_alive=sum(alive[,ns])
    mean_death_period<-mean(dead,na.rm = T)
    d<-data.frame(N1,ns,total_marked,left_alive,mean_death_period)
    results<-results%>%rbind(d)
  }
}

N1<-3000
ns<-20
sim<-ch_sim(N1,ns,coef=MRA_coef)
ch=sim$ch
lengths=sim$lengths
sex=sim$sex
alive<-sim$alive

write.csv(ch,'outputs/sim_ch.csv',row.names = T)
write.csv(lengths,'outputs/lengths_ch.csv',row.names = T)
write.csv(sex,'outputs/sex_ch.csv',row.names = T)

############################################################################
#memoise testing
############################################################################
MRA_coef<-c(3.908292515, -0.006356866,  2.972067636, -0.105286686)
N1=100
ns=15
sim_coef<-MRA_coef#just using coefficients produced by MRA for prior data
iterations<-100
results_mf<-data.frame()
results_obj<-data.frame()
iter_start<-Sys.time()
catch_matrices<-c()
sim<-ch_sim(N1,ns,sim_coef)
catch_hist<-data.frame(sim$ch)
sex_matrix<-sim$sex
lengths_matrix<-sim$lengths

#with memoise
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
  
  sex_orig<-as.matrix(sex_matrix[samp,])
  
  lengths_orig<-as.matrix(lengths_matrix[samp,])
  
  cap_X=lengths_orig
  surv_X=sex_orig
  #but let's also test all the different maximization methods available in optim
  optim_methods<-c("Nelder-Mead", "BFGS")
  iter_df<-data.frame()
  for(i in 1:length(optim_methods)){
    om<-optim_methods[i]
    starttime<-Sys.time()
    CJS_coef<-tryCatch({optim(beta, mf, method = om,hessian=T)$par
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
  }
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
  results_mf<-results_mf%>%rbind(iter_df)
  print(r)
}

#plot results
d<-results_mf
ggplot(d)+
    geom_boxplot(aes(x=om,y=CJS_N))+
    geom_hline(aes(yintercept=N1),color='red')+
    labs(title='with memoise')+
    theme_bw()
ggsave('figures/with_memoise.png')



#without memoise
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
  
  sex_orig<-as.matrix(sex_matrix[samp,])
  
  lengths_orig<-as.matrix(lengths_matrix[samp,])
  
  cap_X=lengths_orig
  surv_X=sex_orig
  #but let's also test all the different maximization methods available in optim
  optim_methods<-c("Nelder-Mead", "BFGS")
  iter_df<-data.frame()
  for(i in 1:length(optim_methods)){
    om<-optim_methods[i]
    starttime<-Sys.time()
    CJS_coef<-tryCatch({optim(beta, CJS_obj, method = om,hessian=T)$par
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
  }
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
  results_obj<-results_obj%>%rbind(iter_df)
  print(r)
}
#plot results
d<-results_obj
ggplot(d)+
  geom_boxplot(aes(x=om,y=CJS_N))+
  geom_hline(aes(yintercept=N1),color='red')+
  labs(title='with memoise')+
  theme_bw()
ggsave('figures/without_memoise.png')

############################################################################
#sim various true N and ns values
############################################################################
MRA_coef<-c(3.908292515, -0.006356866,  2.972067636, -0.105286686)
N1_range<-seq(from=1000,to=1000,by=500)
ns_range<-seq(from=15,to=15,by=5)
library(memoise)
#apparently this memoise package speeds things up, not sure how it works but it does, note that it only speeds up after the initial time the function is called
data<-NULL
schnabel_results<-NULL
for(q in 1:length(N1_range)){
  for(w in 1:length(ns_range)){
    N1=N1_range[q]
    ns=ns_range[w]
    sim_coef<-MRA_coef#just using coefficients produced by MRA for prior data
    iterations<-100
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
      
      lengths_vector<-cap_X[,1]
      sex_vector<-surv_X[,1]
      #but let's also test all the different maximization methods available in optim
      optim_methods<-c('Nelder-Mead',"BFGS","CG")
      iter_df<-data.frame()
      for(i in 3:length(optim_methods)){
        om<-optim_methods[i]
        starttime<-Sys.time()
        CJS_coef<-tryCatch({optim(beta, mf, method = om,hessian=T)$par
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
      }
      #F.cjs_estim
      starttime<-Sys.time()
      #source('scripts/package_testing/F.cjs.estim_testing.R')
      #CJS_coef<-ans$parameters
      #CJS_N<-CJS_escapement(ch,CJS_coef,np)
      #CJS_lik<-ans$loglik
      #endtime<-Sys.time()
      #CJStime<-endtime-starttime
      #d<-data.frame(N1,iter,om='MRA',CJS_N,CJS_lik,CJStime)
      #iter_df<-iter_df%>%rbind(d)
      
      #get schnabel estimate
      #N_schn<-schnabel(ch,ns,nan)[1]
      #CI_schn<-schnabel(ch,ns,nan)[2]
      #schn_data<-data.frame(N_schn,CI_schn)
      
      #iter_df<-iter_df%>%cbind(schn_data)
      results<-results%>%rbind(iter_df)
    }
    results$observations<-ns
    results$total_marked=nrow(ch)
    data<-data%>%rbind(results)
    print(c(paste('N1=',N1),paste('ns=',ns),paste('Total Marked=',nrow(ch))))
    saveRDS(data,'outputs/test_sim_results.rds')
    drop_cache(mf)
  }
}
saveRDS(data,'outputs/test_total_results.rds')
data<-readRDS('outputs/test_total_results.rds')
data$CJS_N<-as.numeric(data$CJS_N)
#confidence intervals
conf.level = 95
alpha = 1 - conf.level/100
lower = alpha/2
upper = 1 - alpha/2
mid=.5
ci<-data
ci<-ci%>%
  group_by(om,N1,observations)%>%
  summarise(lower_ci=ceiling(quantile(CJS_N,probs = c(lower),na.rm=T)),mean=mean(CJS_N,na.rm=T),upper_ci=ceiling(quantile(CJS_N,probs = c(upper),na.rm=T)))

Schn<-unique(select(data,N1,iter,N_schn))
Schn$om<-'Schnabel'
Schn$CJStime<-0
Schn$CJS_lik<-'NA'
Schn<-rename(Schn,CJS_N=N_schn)
data<-select(data,N1,iter,om,CJS_N,CJStime,CJS_lik)
Schn$CJS_lik<-NA
data<-data%>%rbind(Schn)
data$om<-ifelse(data$om=='MRA','Fletcher/MRA',data$om)
data$om<-ifelse(data$om=='Schnabel','Schaefer',data$om)

methods<-unique(ci$om)
#plot results
for(i in 1:length(methods)){
  d<-data%>%filter(om==methods[i])
  ggplot(d)+
    geom_boxplot(aes(x=as.factor(om),y=CJS_N,group=om))+
    facet_wrap(~N1,nrow=1,scales='free_y')+
    geom_hline(aes(yintercept=N1),color='red')+
    labs(title=methods[i])+
    theme_bw()
  ggsave(paste('figures/',methods[i],'_test.png'),width=8,height = 5)
}
#do one plot
times<-data%>%
  group_by(om,N1)%>%
  dplyr::summarise(total_time=sum(CJStime))
N1<-unique(times$N1)

for(i in 1:length(N1)){
  d<-times%>%filter(N1==N1_list[i])
  time1=paste(d$om[1],round(d$total_time[1],2),'min')
  time2=paste(d$om[2],round(d$total_time[2],2),'min')
  time3=paste(d$om[3],round(d$total_time[3],2),'min')
  
  labels[[i]]=paste(time1,time2,time3,sep='\n')
  
}
labels<-unlist(labels)
labels <- as.vector(labels,'character')
times<-data.frame(N1,cat(labels))

ggplot(data)+
  geom_boxplot(aes(x=as.factor(om),y=CJS_N,group=om))+
  facet_wrap(~N1,nrow=1,scales='free_y')+
  geom_text(data=times,inherit.aes = F,aes(x=-Inf,y=Inf,label=labels,hjust='inward',vjust = "inward"),size=2)+
  geom_hline(aes(yintercept=N1),color='red')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('figures/sim_test.png',width=8,height = 4)

#data$CJS_N<-ifelse(data$om=='BFGS' & data$N1==1000,results$CJS_N,data$CJS_N)
