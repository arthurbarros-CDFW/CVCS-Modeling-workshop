#demo script testing
#note: run all chunks in "Modeling_Review.qmd" first
library(escapeMR)
library(dplyr)
library(mra)

###########################################################
#NO CHOPS
############################################################
ans_nochops<-askForData()
CJSscript()

#below is from "CJSScript.R"
chList<-ans_nochops

ch <- ch.for.chat <- chList$ch
sex <- chList$sex
length <- chList$length
covars.exist <- chList$covars.exist
ints <- chList$intervals

#  ---  nj is the total number of carcasses captured (and checked for marks) at sampling occasion j.
nj = apply(ch,2,function(x){ sum(x > 0)})

#  ---  chopped are the number amnong the nj that are chopped (no longer part of the marked or unmarked population)
chopped = apply(ch,2,function(x){sum(x > 1)})

#  ---  rj is the total number of carcasses at occasion j that are released with marks. 
rj = nj - chopped


cap.model = as.formula("~ivar(length,ns)")
surv.model =as.formula("~ivar(sex,ns)")

fit = esc_model(m=7, 
                capt.hist=ch, 
                cap.model=cap.model, 
                surv.model=surv.model, 
                covars.exist = covars.exist,
                sex=sex, 
                length=length, 
                chat=chat, 
                BS=bs, 
                intervals=ints)
#no chops, loglik= -87.14024, N=784

fit_model = F.cjs.estim(cap.model, surv.model, histories=ch, c.hat=chat, intervals= rep(1,ncol(ch)-1))

sex=as.numeric(sex_vector)
length=as.numeric(lengths_vector)
intervals=NULL
if(is.null(intervals)){
  intervals <- rep(1,ncol(ch)-1)
}


#function calls
cap.model = as.formula("~ivar(length,ns)")
surv.model =as.formula("~ivar(sex,ns)")
#cap.model = as.formula("~1")
#surv.model =as.formula("~1")
survival=surv.model
capture=cap.model
histories=as.matrix(ch)
cap.init=NULL
nhat.v.meth=1
sur.init=NULL
c.hat=-1.0
link="logit"
group=NULL
chat=-1
ns = dim(ch)[2]
nan = dim(ch)[1]
Chat_user = chat

control=mra.control()
covars <- F.cr.model.matrix(capture, survival, nan, ns )  
nx <- covars$n.cap.covars  #nx and ny include intercept.  This is total number of parameters
ny <- covars$n.sur.covars

if( length(union( unique(histories), c(0,1,2))) > 3 ) stop("Capture histories must consist of 0's, 1's, and 2's only.")

if( missing(cap.init) ){
  cap.init <- rep(0,nx)
} else if(length(cap.init) < (nx) ){
  cap.init <- c(cap.init, rep(0, nx-length(cap.init)))
} 
if( missing(sur.init) ){
  sur.init <- rep(0,ny)
} else if(length(sur.init) < (ny) ){
  sur.init <- c(sur.init, rep(0, ny-length(sur.init)))
} 
#   Set up the tolerance vector, if not specified, or if not long enough
if( length(control$tol) < (nx+ny) ){
  control$tol <- rep(control$tol, trunc((nx+ny) / length(control$tol))+1)[1:(nx+ny)]
} else if( length(control$tol > (nx+ny)) ){
  control$tol <- control$tol[1:(nx+ny)]
}
if( missing( group )){
  group <- rep(1, nan)
  ng <- 1
} else {
  ng <- length(unique(group))
}
vif <- c.hat

#   Do the estimation, but first allocate room for answers
loglik <- chisq.vif <- df.vif <- 0
parameters <- se.param <- rep(0, nx + ny )
covariance <- matrix( 0, nx+ny, nx+ny )
p.hat <- se.p.hat <- s.hat <- se.s.hat <- matrix( 0, nan, ns )
exit.code <- cov.code <- 0
n.hat <- se.n.hat <- rep(0, ns)
# on entry to .Fortran, maxfn is maximum number of function evals.  On return, maxfn is actual number of evaluations.
maxfn <- control$maxfn  

#   Re-code the link specification to integers
if( link=="logit" ){
  link.code <- 1
} else if( link == "sine" ){
  link.code <- 2
} else if( link == "hazard" ){
  link.code <- 3
} else {
  stop("Unknown link function specified.")
}

if(control$trace) cat( "Calling MRA DLL to maximize likelihood.  Please wait...\n")
ans <- .Fortran( "cjsmod", 
                 nan         = as.integer(nan), 
                 ns          = as.integer(ns), 
                 nx          = as.integer(nx), 
                 ny          = as.integer(ny), 
                 ng          = as.integer(ng), 
                 histories   = as.integer(histories), 
                 group       = as.integer(group), 
                 algorithm   = as.integer(control$algorithm), 
                 cov.meth    = as.integer(control$cov.meth), 
                 link        = as.integer(link.code),
                 nhat.v.meth = as.integer(nhat.v.meth), 
                 capX        = as.double(covars$capX), 
                 survX       = as.double(covars$survX), 
                 cap.init    = as.double(cap.init), 
                 sur.init    = as.double(sur.init), 
                 maxfn       = as.integer(maxfn),
                 beta.tol.vec= as.double(control$tol), 
                 loglik      = as.double(loglik), 
                 vif         = as.double(vif), 
                 chisq.vif   = as.double(chisq.vif), 
                 df.vif      = as.double(df.vif), 
                 parameters  = as.double(parameters),
                 se.param    = as.double(se.param), 
                 covariance  = as.double(covariance), 
                 p.hat       = as.double(p.hat), 
                 se.p.hat    = as.double(se.p.hat), 
                 s.hat       = as.double(s.hat), 
                 se.s.hat    = as.double(se.s.hat), 
                 n.hat       = as.double(n.hat), 
                 se.n.hat    = as.double(se.n.hat), 
                 exit.code   = as.integer(exit.code), 
                 cov.code    = as.integer(cov.code), 
                 intervals   = as.double(intervals), 
                 PACKAGE="mra" 
) 


#COMPARE TO MY CODE
ch<-read.csv('data/demo_CH.csv')
covars<-read.csv('data/demo_Covariates.csv')

ch=ch[-1] #remove disctag vector from capture histories
covars$sex<-ifelse(covars$sex=='F',1,0) #change sex to numeric value 

sex_matrix<-matrix(covars$sex,nrow=nrow(ch),ncol=ncol(ch))
lengths_matrix<-matrix(covars$length,nrow=nrow(ch),ncol=ncol(ch))

mra_beta<-ans$parameters
#3.012170268 -0.007657378  1.343169594  0.748672276

CJS_loglik(beta=mra_beta,cap_X=lengths_matrix,surv_X=sex_matrix,ch=ch)
#no chops, loglik= -87.14024 which is good



###########################################################
#WITH CHOPS
############################################################
ans_withchops<-askForData()

#below is from "CJSScript.R"
chList<-ans_withchops

ch <- ch.for.chat <- chList$ch
sex <- chList$sex
length <- chList$length
covars.exist <- chList$covars.exist
ints <- chList$intervals

#  ---  nj is the total number of carcasses captured (and checked for marks) at sampling occasion j.
nj = apply(ch,2,function(x){ sum(x > 0)})

#  ---  chopped are the number amnong the nj that are chopped (no longer part of the marked or unmarked population)
chopped = apply(ch,2,function(x){sum(x > 1)})

#  ---  rj is the total number of carcasses at occasion j that are released with marks. 
rj = nj - chopped


cap.model = as.formula("~ivar(length,ns)")
surv.model =as.formula("~ivar(sex,ns)")

estim.chat = escapeMR::esc_model(capt.hist=ch.for.chat, 
                                 cap.model=as.formula("~ 1"), 
                                 surv.model=as.formula("~ 1"), 
                                 intervals=ints, 
                                 chat=-1, 
                                 full.model=TRUE, 
                                 BS=FALSE)
chat = ifelse(estim.chat >= 1.1, estim.chat, 1)
bs = FALSE
fm = FALSE
ns = dim(ch)[2]
nan = dim(ch)[1]

fit = esc_model(m=7, 
                capt.hist=ch, 
                cap.model=cap.model, 
                surv.model=surv.model, 
                covars.exist = covars.exist,
                sex=sex, 
                length=length, 
                chat=chat, 
                BS=bs, 
                intervals=ints)
#with chops, loglik= -87.14024, the same as without chops? N= 2116

fit_model = F.cjs.estim(cap.model, surv.model, histories=ch, c.hat=chat, intervals= rep(1,ncol(ch)-1))

sex=as.numeric(sex)
length=as.numeric(sex)
intervals=NULL
if(is.null(intervals)){
  intervals <- rep(1,ncol(ch)-1)
}


#function calls
cap.model = as.formula("~ivar(length,ns)")
surv.model =as.formula("~ivar(sex,ns)")
#cap.model = as.formula("~1")
#surv.model =as.formula("~1")
survival=surv.model
capture=cap.model
histories=as.matrix(ch)
cap.init=NULL
nhat.v.meth=1
sur.init=NULL
c.hat=-1.0
link="logit"
group=NULL
chat=-1
ns = dim(ch)[2]
nan = dim(ch)[1]
Chat_user = chat

control=mra.control()

covars <- F.cr.model.matrix(capture, survival, nan, ns )  


nx <- covars$n.cap.covars  #nx and ny include intercept.  This is total number of parameters
ny <- covars$n.sur.covars

if( length(union( unique(histories), c(0,1,2))) > 3 ) stop("Capture histories must consist of 0's, 1's, and 2's only.")

if( missing(cap.init) ){
  cap.init <- rep(0,nx)
} else if(length(cap.init) < (nx) ){
  cap.init <- c(cap.init, rep(0, nx-length(cap.init)))
} 
if( missing(sur.init) ){
  sur.init <- rep(0,ny)
} else if(length(sur.init) < (ny) ){
  sur.init <- c(sur.init, rep(0, ny-length(sur.init)))
} 
#   Set up the tolerance vector, if not specified, or if not long enough
if( length(control$tol) < (nx+ny) ){
  control$tol <- rep(control$tol, trunc((nx+ny) / length(control$tol))+1)[1:(nx+ny)]
} else if( length(control$tol > (nx+ny)) ){
  control$tol <- control$tol[1:(nx+ny)]
}
if( missing( group )){
  group <- rep(1, nan)
  ng <- 1
} else {
  ng <- length(unique(group))
}
vif <- c.hat

#   Do the estimation, but first allocate room for answers
loglik <- chisq.vif <- df.vif <- 0
parameters <- se.param <- rep(0, nx + ny )
covariance <- matrix( 0, nx+ny, nx+ny )
p.hat <- se.p.hat <- s.hat <- se.s.hat <- matrix( 0, nan, ns )
exit.code <- cov.code <- 0
n.hat <- se.n.hat <- rep(0, ns)
# on entry to .Fortran, maxfn is maximum number of function evals.  On return, maxfn is actual number of evaluations.
maxfn <- control$maxfn  

#   Re-code the link specification to integers
if( link=="logit" ){
  link.code <- 1
} else if( link == "sine" ){
  link.code <- 2
} else if( link == "hazard" ){
  link.code <- 3
} else {
  stop("Unknown link function specified.")
}

if(control$trace) cat( "Calling MRA DLL to maximize likelihood.  Please wait...\n")
ans <- .Fortran( "cjsmod", 
                 nan         = as.integer(nan), 
                 ns          = as.integer(ns), 
                 nx          = as.integer(nx), 
                 ny          = as.integer(ny), 
                 ng          = as.integer(ng), 
                 histories   = as.integer(histories), 
                 group       = as.integer(group), 
                 algorithm   = as.integer(control$algorithm), 
                 cov.meth    = as.integer(control$cov.meth), 
                 link        = as.integer(link.code),
                 nhat.v.meth = as.integer(nhat.v.meth), 
                 capX        = as.double(covars$capX), 
                 survX       = as.double(covars$survX), 
                 cap.init    = as.double(cap.init), 
                 sur.init    = as.double(sur.init), 
                 maxfn       = as.integer(maxfn),
                 beta.tol.vec= as.double(control$tol), 
                 loglik      = as.double(loglik), 
                 vif         = as.double(vif), 
                 chisq.vif   = as.double(chisq.vif), 
                 df.vif      = as.double(df.vif), 
                 parameters  = as.double(parameters),
                 se.param    = as.double(se.param), 
                 covariance  = as.double(covariance), 
                 p.hat       = as.double(p.hat), 
                 se.p.hat    = as.double(se.p.hat), 
                 s.hat       = as.double(s.hat), 
                 se.s.hat    = as.double(se.s.hat), 
                 n.hat       = as.double(n.hat), 
                 se.n.hat    = as.double(se.n.hat), 
                 exit.code   = as.integer(exit.code), 
                 cov.code    = as.integer(cov.code), 
                 intervals   = as.double(intervals), 
                 PACKAGE="mra" 
) 

#loglik= -87.54393

#COMPARE TO MY CODE
ch<-read.csv('data/demo_CH.csv')
covars<-read.csv('data/demo_Covariates.csv')
chops<-read.csv('data/demo_Chops.csv')

prepped_data<-CJS_data_prep(ch,chops,covars)

ch_prepped<-prepped_data$ch
cap_X_prepped<-prepped_data$lengths_matrix
surv_X_prepped<-prepped_data$sex_matrix

mra_beta<-ans$parameters
#-2.9248039  0.9582045  2.3739011 -0.8700545

CJS_loglik(beta=ans$parameters,cap_X=cap_X_prepped,surv_X=surv_X_prepped,ch=ch_prepped)
#with chops, getting NAN?

#DETERMINE THAT BECAUSE ASSIGNMENT of covars to chop data is random, difficult to replicate exact outputs

