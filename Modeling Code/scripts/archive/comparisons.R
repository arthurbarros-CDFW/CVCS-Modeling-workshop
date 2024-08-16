#model comparisons
rm( list = ls()) #clear env
library(dplyr)
library(stats4)
library(escapeMR)
library(ggplot2)
library(mra)

source('scripts/CJS_functions.R')
source('scripts/schnabel_est.R')

############################################################################
#First we compare with a real dataset
############################################################################

#read data for my script
ch<-read.csv('data/demo_CH.csv')
ch<-ch[-1]
covars<-read.csv('data/demo_Covariates.csv')
sex<-covars$Sex
sex<-ifelse(sex=='M',1,0)
lengths<-covars$Length

#function calls
cap.model = as.formula("~1")
surv.model =as.formula("~ivar(length,ns)")

#Set Global Variables, can't figure out why I can't pass these in the function...
nan=nrow(ch)
ns=ncol(ch)
p=np=4
nx=ny=p/2
beta= numeric(p)

lengths_matrix<-matrix(lengths)
lengths_matrix<-lengths_matrix[, rep(1,each=ns)]

sex_matrix<-matrix(sex)
sex_matrix<-sex_matrix[, rep(1,each=ns)]
cap_X=lengths_matrix
surv_X=sex_matrix

#define tolerances
beta_tol<-c(1e-07,1e-07,1e-07,1e-07)
central_diffs=TRUE
mode=1

CJS_coef<-optim(beta, CJS_obj, method = "BFGS",hessian=T)$par
CJS_N<-CJS_escapement(ch,CJS_coef,np) #24.20663
CJS_lik<-CJS_loglik(CJS_coef,p) #-13.82324, slightly lower than what MRA produces
endtime<-Sys.time()


#F.cjs_estim
starttime<-Sys.time()
source('scripts/F.cjs.estim_testing.R')
MRA_coef<-ans$parameters
MRA_N<-CJS_escapement(ch,MRA_coef,np)#24.82458
endtime<-Sys.time()
MRAtime<-endtime-starttime #~0.05 sec

#Let's just prove that log likelihood is the same when we set the coef for my CJS estimates to what MRA produces
CJS_lik<-CJS_loglik(MRA_coef,p) #-13.82954
MRA_lik<-ans$loglik #-13.82954

#Let's double check our results using the escapeMR shell script
CJSscript() #Messages:
#model: 7:   capture related to length and survival related to sex
#Log likelihood =  -13.8295429418454
#reported escapement=25
escapeMR()

############################################################################
#bootstrapping FR 2020
############################################################################
#function calls
cap.model = as.formula("~1")
surv.model =as.formula("~ivar(length,ns)")
#read data for my script
catch_hist<-read.csv('data/FR20-9_3_Final Capture History Matrix.csv')
chops<-read.csv('data/FR20-9_4_Final Chops Matrix Table.csv')
covariates<-read.csv('data/FR20-9_5_Final Covariate Table.csv')
covariates[,2]<-ifelse(covariates[,2]=='FALSE','F',covariates[,2])
iterations<-100
results<-data.frame()

#prep in chops data
catch_hist<-catch_hist[-1]
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
colnames(clean_chops)<-colnames(catch_hist)
#add in chops
catch_hist<-catch_hist%>%rbind(clean_chops)

catch_matrices<-c()
library(memoise)
mf<-memoise(CJS_obj) #apparently this memoise package speeds things up, not sure how it works but it does, note that it only speeds up after the initial time the function is called
iter_start<-Sys.time()
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
  
  #if there is chopped data, we will end up with some records not having corresponding covariates
  #we can generate covariates for those chopped histories based on averages
  set.seed(83842912)  # so that random assignments of sex and length is repeatable
  covariates_orig<-as.matrix(covariates[samp,])
  
  sex_vector<-covariates_orig[,2]
  sex_vector<-ifelse(sex_vector=='F',1,0)
  no.miss = sex_vector[!is.na(sex_vector)]
  prop.female<-mean(no.miss)
  prop.male<-1-prop.female
  sex_vector<-ifelse(is.na(sex_vector),sample(c(1,0),1,prob=c(prop.female,prop.male)),sex_vector)
  sex_matrix<-matrix(sex_vector,nrow=nan,ncol=ns)
  
  lengths_vector<-covariates_orig[,3]
  females<-covariates_orig[covariates_orig[,2]=='F',]
  males<-covariates_orig[covariates_orig[,2]=='M',]
  avg.female.length<-round(mean(as.numeric(females[,3]),na.rm=T),1)
  avg.male.length<-round(mean(as.numeric(males[,3]),na.rm=T),1)
  for(i in 1:length(lengths_vector)){
    if(is.na(lengths_vector[i])){
      lengths_vector[i] = ifelse(sex_matrix[i,1] == 1,  avg.female.length, avg.male.length)
    }
  }
  lengths_matrix<-matrix(as.numeric(lengths_vector),nrow=nan,ncol=ns)
  
  cap_X=lengths_matrix
  surv_X=lengths_matrix
  
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
    d<-data.frame(iter,om,CJS_N,CJS_lik,CJStime)
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
  d<-data.frame(iter,om='MRA',CJS_N,CJS_lik,CJStime)
  iter_df<-iter_df%>%rbind(d)
  results<-results%>%rbind(iter_df)
  catch_matrices[[r]]<-ch
  print(r)
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
  geom_text(data=ci,aes(x=mid_ci,y=17,label=paste('Mean:',mid_ci)))+
  geom_text(data=ci,aes(x=lower_ci,y=17,label=paste('Lower:',lower_ci)))+
  geom_text(data=ci,aes(x=upper_ci,y=17,label=paste('Upper:',upper_ci)))+
  theme_classic()
ggsave('figures/FR2020_bootstrap.png',width=8,height = 5)
saveRDS(data,'outputs/WR2020.rds')

############################################################################
#bootstrapping WR 2020
############################################################################

#read data for my script
catch_hist<-read.csv('data/WR20-9_3_Final Capture History Matrix.csv')
chops<-read.csv('data/WR20-9_4_Final Chops Matrix Table.csv')
covariates<-read.csv('data/WR20-9_5_Final Covariate Table.csv')
covariates[,2]<-ifelse(covariates[,2]=='FALSE','F',covariates[,2])
iterations<-1
results<-data.frame()

#prep in chops data
catch_hist<-catch_hist[-1]
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
colnames(clean_chops)<-colnames(catch_hist)
#add in chops
catch_hist<-catch_hist%>%rbind(clean_chops)

catch_matrices<-c()
library(memoise)
mf<-memoise(CJS_obj) #apparently this memoise package speeds things up, not sure how it works but it does, note that it only speeds up after the initial time the function is called
iter_start<-Sys.time()
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
  
  #if there is chopped data, we will end up with some records not having corresponding covariates
  #we can generate covariates for those chopped histories based on averages
  set.seed(83842912)  # so that random assignments of sex and length is repeatable
  covariates_orig<-as.matrix(covariates[samp,])
  
  sex_vector<-covariates_orig[,2]
  sex_vector<-ifelse(sex_vector=='F',1,0)
  no.miss = sex_vector[!is.na(sex_vector)]
  prop.female<-mean(no.miss)
  prop.male<-1-prop.female
  sex_vector<-ifelse(is.na(sex_vector),sample(c(1,0),1,prob=c(prop.female,prop.male)),sex_vector)
  sex_matrix<-matrix(sex_vector,nrow=nan,ncol=ns)
  
  lengths_vector<-covariates_orig[,3]
  females<-covariates_orig[covariates_orig[,2]=='F',]
  males<-covariates_orig[covariates_orig[,2]=='M',]
  avg.female.length<-round(mean(as.numeric(females[,3]),na.rm=T),1)
  avg.male.length<-round(mean(as.numeric(males[,3]),na.rm=T),1)
  for(i in 1:length(lengths_vector)){
    if(is.na(lengths_vector[i])){
      lengths_vector[i] = ifelse(sex_matrix[i,1] == 1,  avg.female.length, avg.male.length)
    }
  }
  lengths_matrix<-matrix(as.numeric(lengths_vector),nrow=nan,ncol=ns)
  
  cap_X=lengths_matrix
  surv_X=sex_matrix
  
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
    d<-data.frame(iter,om,CJS_N,CJS_lik,CJStime)
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
  d<-data.frame(iter,om='MRA',CJS_N,CJS_lik,CJStime)
  iter_df<-iter_df%>%rbind(d)
  results<-results%>%rbind(iter_df)
  catch_matrices[[r]]<-ch
  print(r)
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
saveRDS(data,'outputs/WR2020.rds')
data<-readRDS('outputs/WR2020.rds')
ggplot(data,aes(x=CJS_N))+
  geom_histogram(color = "#000000", fill = "#0099F8")+
  facet_grid(factor(om,levels=c('Nelder-Mead','BFGS','MRA'))~.)+
  geom_segment(data=ci,aes(x=lower_ci,xend=lower_ci,y=0,yend=60),linewidth=1,linetype='dashed')+
  geom_segment(data=ci,aes(x=upper_ci,xend=upper_ci,y=0,yend=60),linewidth=1,linetype='dashed')+
  geom_segment(data=ci,aes(x=mid_ci,xend=mid_ci,y=0,yend=60),linewidth=1,linetype='dashed')+
  #geom_text(data=ci,aes(x=mid_ci,y=17,label=paste('Mean:',mid_ci)))+
  #geom_text(data=ci,aes(x=lower_ci,y=17,label=paste('Lower:',lower_ci)))+
  #geom_text(data=ci,aes(x=upper_ci,y=17,label=paste('Upper:',upper_ci)))+
  theme_classic()
ggsave('figures/WR2020_bootstrap.png',width=8,height = 5)


############################################################################
#bootstrapping FR 2021
############################################################################

#read data for my script
catch_hist<-read.csv('data/fall-21 9_3_Final Capture History Matrix.csv')
chops<-read.csv('data/fall-21 9_4_Final Chops Matrix Table.csv')
covariates<-read.csv('data/fall-21 9_5_Final Covariate Table.csv')
iterations<-1
results<-data.frame()

#prep in chops data
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
colnames(clean_chops)<-colnames(catch_hist[-1])

catch_matrices<-c()
library(memoise)
mf<-memoise(CJS_obj) #apparently this memoise package speeds things up, not sure how it works but it does, note that it only speeds up after the initial time the function is called
iter_start<-Sys.time()
for(r in 1:iterations){
  iter<-r
  index = 1:dim(catch_hist)[1]
  samp<- sample(index, replace = T)
  ch = catch_hist[samp,]
  ch<-ch[-1]
  #Set Global Variables, can't figure out why I can't pass these in the function...
  #add in chops
  ch<-ch%>%rbind(clean_chops)
  nan=nrow(ch)
  ns=ncol(ch)
  p=np=4
  nx=ny=p/2
  beta= numeric(p)
  
  sex_orig<-as.matrix(covariates[samp,])
  sex_matrix<-sex_orig[,2]
  sex_matrix<-matrix(sex_matrix,nrow=nan,ncol=ns)
  
  lengths_orig<-as.matrix(covariates[samp,])
  lengths_matrix<-lengths_orig[,3]
  lengths_matrix<-matrix(lengths_matrix,nrow=nan,ncol=ns)
  
  cap_X=lengths_matrix
  surv_X=sex_matrix
  
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
    d<-data.frame(iter,om,CJS_N,CJS_lik,CJStime)
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
  d<-data.frame(iter,om='MRA',CJS_N,CJS_lik,CJStime)
  iter_df<-iter_df%>%rbind(d)
  results<-results%>%rbind(iter_df)
  catch_matrices[[r]]<-ch
  print(r)
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

#plot results
data<-results
data$CJS_N<-as.numeric(data$CJS_N)
saveRDS(data,'outputs/FR2021.rds')
data<-readRDS('outputs/FR2021.rds')
ggplot(data,aes(x=CJS_N))+
  geom_histogram(color = "#000000", fill = "#0099F8")+
  facet_grid(factor(om,levels=c('Nelder-Mead','BFGS','MRA'))~.)+
  #geom_segment(data=ci,aes(x=lower_ci,xend=lower_ci,y=0,yend=15),linewidth=1,linetype='dashed')+
  #geom_segment(data=ci,aes(x=upper_ci,xend=upper_ci,y=0,yend=15),linewidth=1,linetype='dashed')+
  #geom_segment(data=ci,aes(x=mid_ci,xend=mid_ci,y=0,yend=15),linewidth=1,linetype='dashed')+
  #geom_text(data=ci,aes(x=mid_ci,y=17,label=paste('Mean:',mid_ci)))+
  #geom_text(data=ci,aes(x=lower_ci,y=17,label=paste('Lower:',lower_ci)))+
  #geom_text(data=ci,aes(x=upper_ci,y=17,label=paste('Upper:',upper_ci)))+
  theme_classic()
ggsave('figures/FR2021_bootstrap.png',width=8,height = 5)

############################################################################
#bootstrapping WR 2021
############################################################################

#read data for my script
catch_hist<-read.csv('data/9_3_Final Capture History Matrix.csv')
chops<-read.csv('data/9_4_Final Chops Matrix Table.csv')
covariates<-read.csv('data/9_5_Final Covariate Table.csv')
iterations<-100
results<-data.frame()

#prep in chops data
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
colnames(clean_chops)<-colnames(catch_hist[-1])

catch_matrices<-c()
library(memoise)
mf<-memoise(CJS_obj) #apparently this memoise package speeds things up, not sure how it works but it does, note that it only speeds up after the initial time the function is called
iter_start<-Sys.time()
for(r in 1:iterations){
  iter<-r
  index = 1:dim(catch_hist)[1]
  samp<- sample(index, replace = T)
  ch = catch_hist[samp,]
  ch<-ch[-1]
  #Set Global Variables, can't figure out why I can't pass these in the function...
  #add in chops
  ch<-ch%>%rbind(clean_chops)
  nan=nrow(ch)
  ns=ncol(ch)
  p=np=4
  nx=ny=p/2
  beta= numeric(p)
  
  sex_orig<-as.matrix(covariates[samp,])
  sex_matrix<-sex_orig[,2]
  sex_matrix<-matrix(sex_matrix,nrow=nan,ncol=ns)
  
  lengths_orig<-as.matrix(covariates[samp,])
  lengths_matrix<-lengths_orig[,3]
  lengths_matrix<-matrix(lengths_matrix,nrow=nan,ncol=ns)
  
  cap_X=lengths_matrix
  surv_X=sex_matrix
  
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
    d<-data.frame(iter,om,CJS_N,CJS_lik,CJStime)
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
  d<-data.frame(iter,om='MRA',CJS_N,CJS_lik,CJStime)
  iter_df<-iter_df%>%rbind(d)
  results<-results%>%rbind(iter_df)
  catch_matrices[[r]]<-ch
  print(r)
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

#plot results
data<-results
data$CJS_N<-as.numeric(data$CJS_N)
ggplot(data,aes(x=CJS_N))+
  geom_histogram(color = "#000000", fill = "#0099F8")+
  facet_grid(factor(om,levels=c('Nelder-Mead','BFGS','MRA'))~.)+
  geom_segment(data=ci,aes(x=lower_ci,xend=lower_ci,y=0,yend=15),linewidth=1,linetype='dashed')+
  geom_segment(data=ci,aes(x=upper_ci,xend=upper_ci,y=0,yend=15),linewidth=1,linetype='dashed')+
  geom_segment(data=ci,aes(x=mid_ci,xend=mid_ci,y=0,yend=15),linewidth=1,linetype='dashed')+
  geom_text(data=ci,aes(x=mid_ci,y=17,label=paste('Mean:',mid_ci)))+
  geom_text(data=ci,aes(x=lower_ci,y=17,label=paste('Lower:',lower_ci)))+
  geom_text(data=ci,aes(x=upper_ci,y=17,label=paste('Upper:',upper_ci)))+
  theme_classic()
saveRDS(data,'outputs/WR2021.rds')
ggsave('figures/WR2021_bootstrap.png',width=8,height = 5)
