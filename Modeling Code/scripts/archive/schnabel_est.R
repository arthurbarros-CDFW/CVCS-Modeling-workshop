############################################################################
#establish Schnabel estimate function
############################################################################

schnabel<-function(ch,ns,nan){
  source('scripts/package_testing/CJS_functions.R')
  C<-NA
  M<-NA
  chops<-NA
  R<-NA
  U<-NA
  Mt<-NA
  first<-location(nan,ns,ch)$first
  last<-location(nan,ns,ch)$last
  recaps<-matrix(NA,nan,ns)
  for(i in 1:nan){
    for(j in 1:ns){
      recaps[i,j]<-ifelse(j>first[i] & ch[i,j]>=1,1,0)
    }
  }
  recaps<-as.data.frame(recaps)
  for(j in 1:ns){
    C[j]<-nrow(ch[ch[j]>=1,])
    M[j]<-nrow(ch[ch[j]==1,])
    chops[j]<-nrow(ch[ch[j]==2,])
    R[j]<-nrow(recaps[recaps[j]==1,])
  }
  for(j in 2:ns){
    Mt[j]<-sum(M[2:j])
  }
  Mt[1]=0
  Nest<-sum(Mt*C)/(sum(R))
  Nvar<-sum(R)/(sum(Mt*C)^2)
  SE<-sqrt(Nvar)
  
  t_crit=qt(p=0.05,df=(ns-1),lower.tail = F)
  margin.error<-t_crit*SE
  CI<-1/(t_crit*SE)
  return(c(Nest,CI))
}


