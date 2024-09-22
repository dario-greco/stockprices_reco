#
# Code for reproducing results in Section 6
#
###########################################################

rm(list=ls())

library(rio)
library(here)
library(sandwich)
library(lmtest)
library(rugarch)

setwd(here("RDSForecast"))

SR.test <- function(Price, predPrice, pred2) {
  
  ret <- NULL
  for (i in 1:length(Price)) {
    ret[i] <- (Price[i]-Price[i-1])/Price[i-1]
  }
  ret <- na.omit(ret)
  
  rethat <- NULL
  for (i in 1:length(Price)) {
    rethat[i] <- (predPrice[i]-Price[i-1])/Price[i-1]
  }
  rethat <- na.omit(rethat)
    
  rethat2 <- NULL
  for (i in 1:length(Price)) {
      rethat2[i] <- (pred2[i]-Price[i-1])/Price[i-1]
    }
  rethat2 <- na.omit(rethat2)
    
  signx = sign(ret)
  sign1 = sign(rethat)
  sign1[sign1 == 0] = 1
  sign2 = sign(rethat2)
  sign2[sign2 == 0] = 1
    
  r_t = sign1 * ret
  r_t2 = sign2 * ret
    
  stat = ((mean(r_t)/sd(r_t)) - (mean(r_t2)/sd(r_t2))) * 100
    
  Shr = sd(r_t2)*mean(r_t) - sd(r_t)*mean(r_t2)
  Theta = (var(r_t)*var(r_t2)/length(ret)) * ( 2*cov(r_t,r_t2) + 0.5 * ( (mean(r_t)/sd(r_t))^2 + (mean(r_t2)/sd(r_t2))^2 - ((mean(r_t)/sd(r_t))*(mean(r_t2)/sd(r_t2)) * (1+cov(r_t, r_t2)^2)) ) )
    
  p.val = 1-pnorm(Shr/sqrt(Theta))
    
  return(list=c("Statistic"=stat, p.val=p.val))
  
}

Ytotlist <- import("Ytotlist.RDS")

Ytot.I <- Ytotlist[[1]]
Ytot.IND <- Ytotlist[[2]]
Ytot.MRKT <- Ytotlist[[3]]
Ytot.EUCL <- Ytotlist[[4]]
Ytot.COR <- Ytotlist[[5]]
Ytot.ARMA <- Ytotlist[[6]]
Ytot.all <- Ytotlist[[7]]

train.l <- 400
out.l <- nrow(Ytot.all)-train.l
h <- c(1,3,6,12) 
type <- c("for.h","MinTrec.h")

InvRes <- list()
InvPval <- list()

for (l in 1:4){
  
  print(paste0("Computing",l,"/4..."))
  
  IND.f<-import(paste0(type[1],h[l],"IND.RDS"))
  IND.rec<-import(paste0(type[2],h[l],"IND.RDS"))
  MRKT.f<-import(paste0(type[1],h[l],"MRKT.RDS"))
  MRKT.rec<-import(paste0(type[2],h[l],"MRKT.RDS"))
  EUCL.f<-import(paste0(type[1],h[l],"EUCL.RDS"))
  EUCL.rec<-import(paste0(type[2],h[l],"EUCL.RDS"))
  COR.f<-import(paste0(type[1],h[l],"COR.RDS"))
  COR.rec<-import(paste0(type[2],h[l],"COR.RDS"))
  ARMA.f<-import(paste0(type[1],h[l],"ARMA.RDS"))
  ARMA.rec<-import(paste0(type[2],h[l],"ARMA.RDS"))
  ALL.f<-import(paste0(type[1],h[l],"ALL.RDS"))
  ALL.rec<-import(paste0(type[2],h[l],"ALL.RDS"))
  NOCLUST.f<-import(paste0(type[1],h[l],"NOCLUST.RDS"))
  NOCLUST.rec<-import(paste0(type[2],h[l],"NOCLUST.RDS"))
  BU.rec<-import(paste0("rec.h",h[l],"BU",".RDS"))
  RW.f <- import(paste0("forRW.h",h[l],".RDS"))
  
  yval<-vector("numeric",length(IND.f))
  for (i in 1:length(yval)) {
    yval[i]<-Ytot.all[(train.l+i+h[l]-1),1]
  }
  
  yhat.IND <- NULL
  for (i in 1:length(IND.f)) {
    yhat.IND[i] <- sum(IND.rec[[i]][h[l],])
  }
  
  yhat.MRKT <- NULL
  for (i in 1:length(IND.f)) {
    yhat.MRKT[i] <- sum(MRKT.rec[[i]][h[l],])
  }
  
  yhat.EUCL <- NULL
  for (i in 1:length(IND.f)) {
    yhat.EUCL[i] <- sum(EUCL.rec[[i]][h[l],])
  }
  
  yhat.COR <- NULL
  for (i in 1:length(IND.f)) {
    yhat.COR[i] <- sum(COR.rec[[i]][h[l],])
  }
  
  yhat.ARMA <- NULL
  for (i in 1:length(IND.f)) {
    yhat.ARMA[i] <- sum(ARMA.rec[[i]][h[l],])
  }
  
  yhat.ALL <- NULL
  for (i in 1:length(IND.f)) {
    yhat.ALL[i] <- ifelse(h[l]==1, sum(ALL.rec[[i]]), sum(ALL.rec[[i]][,h[l]]))
  }
  
  yhat.BU <- NULL
  for (i in 1:length(IND.f)) {
    yhat.BU[i] <- sum(BU.rec[[i]][-1,h[l]])
  }
  
  yhat.NOCLUST <- NULL
  for (i in 1:length(IND.f)) {
    yhat.NOCLUST[i] <- sum(NOCLUST.rec[[i]][h[l],])
  }
  
  yhat.Base <- NULL
  for (i in 1:length(IND.f)) {
    yhat.Base[i] <- NOCLUST.f[[i]][[1]][h[l],1]
  }
  
  
  # Proportion test: ARMA vs reconciliation h=1
  
  InvMat1 <- matrix(NA,6,3)
  InvMat2 <- matrix(NA,6,3)
  
  InvMat1[1,1] <- SR.test(yval, yhat.IND, yhat.Base)[1]
  InvMat2[1,1] <- SR.test(yval, yhat.IND, yhat.Base)[2]
  InvMat1[2,1] <- SR.test(yval, yhat.MRKT, yhat.Base)[1]
  InvMat2[2,1] <- SR.test(yval, yhat.MRKT, yhat.Base)[2]
  InvMat1[3,1] <- SR.test(yval, yhat.EUCL, yhat.Base)[1]
  InvMat2[3,1] <- SR.test(yval, yhat.EUCL, yhat.Base)[2]
  InvMat1[4,1] <- SR.test(yval, yhat.COR, yhat.Base)[1]
  InvMat2[4,1] <- SR.test(yval, yhat.COR, yhat.Base)[2]
  InvMat1[5,1] <- SR.test(yval, yhat.ARMA, yhat.Base)[1]
  InvMat2[5,1] <- SR.test(yval, yhat.ARMA, yhat.Base)[2]
  InvMat1[6,1] <- SR.test(yval, yhat.ALL, yhat.Base)[1]
  InvMat2[6,1] <- SR.test(yval, yhat.ALL, yhat.Base)[2]
  
  InvMat1[1,2] <- SR.test(yval, yhat.IND, yhat.BU)[1]
  InvMat2[1,2] <- SR.test(yval, yhat.IND, yhat.BU)[2]
  InvMat1[2,2] <- SR.test(yval, yhat.MRKT, yhat.BU)[1]
  InvMat2[2,2] <- SR.test(yval, yhat.MRKT, yhat.BU)[2]
  InvMat1[3,2] <- SR.test(yval, yhat.EUCL, yhat.BU)[1]
  InvMat2[3,2] <- SR.test(yval, yhat.EUCL, yhat.BU)[2]
  InvMat1[4,2] <- SR.test(yval, yhat.COR, yhat.BU)[1]
  InvMat2[4,2] <- SR.test(yval, yhat.COR, yhat.BU)[2]
  InvMat1[5,2] <- SR.test(yval, yhat.ARMA, yhat.BU)[1]
  InvMat2[5,2] <- SR.test(yval, yhat.ARMA, yhat.BU)[2]
  InvMat1[6,2] <- SR.test(yval, yhat.ALL, yhat.BU)[1]
  InvMat2[6,2] <- SR.test(yval, yhat.ALL, yhat.BU)[2]
  
  InvMat1[1,3] <- SR.test(yval, yhat.IND, yhat.NOCLUST)[1]
  InvMat2[1,3] <- SR.test(yval, yhat.IND, yhat.NOCLUST)[2]
  InvMat1[2,3] <- SR.test(yval, yhat.MRKT, yhat.NOCLUST)[1]
  InvMat2[2,3] <- SR.test(yval, yhat.MRKT, yhat.NOCLUST)[2]
  InvMat1[3,3] <- SR.test(yval, yhat.EUCL, yhat.NOCLUST)[1]
  InvMat2[3,3] <- SR.test(yval, yhat.EUCL, yhat.NOCLUST)[2]
  InvMat1[4,3] <- SR.test(yval, yhat.COR, yhat.NOCLUST)[1]
  InvMat2[4,3] <- SR.test(yval, yhat.COR, yhat.NOCLUST)[2]
  InvMat1[5,3] <- SR.test(yval, yhat.ARMA, yhat.NOCLUST)[1]
  InvMat2[5,3] <- SR.test(yval, yhat.ARMA, yhat.NOCLUST)[2]
  InvMat1[6,3] <- SR.test(yval, yhat.ALL, yhat.NOCLUST)[1]
  InvMat2[6,3] <- SR.test(yval, yhat.ALL, yhat.NOCLUST)[2]
  
  InvRes[[l]] <- InvMat1
  InvPval[[l]] <- InvMat2
  
}

print(xtable::xtable(InvRes[[3]], digits=2), include.rownames = F)
print(xtable::xtable(InvRes[[4]], digits=2), include.rownames = F)
round(InvPval[[4]],2)