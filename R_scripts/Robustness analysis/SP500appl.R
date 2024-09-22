#
# Code for reproducing results in Section 6 (SP500)
#
###########################################################

rm(list=ls())

library(forecast)
library(hts)
library(TSclust)
library(astsa)
library(BatchGetSymbols)
library(rio)
library(cluster)
library(foreach)
library(sandwich)
library(lmtest)
library(doParallel)
library(ggcorrplot)
library(here)
library(reshape2)

# Load data:

SP.dat <- rio::import("sp500_companies.csv")# Tickers, exchange, industry
DataP0_wide <- rio::import("sp500wide.RData") # Prices SP500
Basey <- DataP0_wide[,-c(which(apply(DataP0_wide, 2, function(x){sum(is.na(x))})!=0))]
Basey <- Basey[,-1]
SP.dat <- SP.dat[SP.dat$Symbol %in% colnames(Basey), ]

# Plotting

DJI_df <- DataP0_wide |>
  mutate(
    DJI = rowSums(Basey)
  ) |>
  #select(-Date) |>
  tidyr::pivot_longer(-Date, names_to = "Stock", values_to = "Price")


stocks <- sort(unique(DJI_df$Stock))
stocks <- c(stocks[stocks != "DJI"], "DJI")
DJI_df |>
  mutate(Index = if_else(Stock == "DJI", "S&P500", "Stocks")) |>
  ggplot(aes(x = Date, y = log(Price), group = Stock, col = Stock)) +
  geom_line() +
  facet_grid(Index ~ ., scales = "free_y") +
  labs(title = "S&P500 Index with Component Stocks") +
  scale_x_date(date_breaks = "1 year", minor_breaks = "1 month", date_labels = "%Y") +
  theme_bw() +
  labs(y = "Price") +
  scale_color_manual(breaks = stocks, values = c(scales::hue_pal()(length(stocks)-1),"#000000")) +
  guides(col = "none")

# A. Discovering clustering structure:

# 0. No clustering:

Hy.I <- hts(Basey)
Ytot.I <- aggts(Hy.I)

# 1. Subjective clusters:

# Industry-based clusters:

nclust.IND <- length(unique(SP.dat$Sector))
SP.dat2 <- SP.dat[order(SP.dat$Symbol),]
clust.IND <- as.numeric(as.factor(SP.dat$Sector))
SP.dat2 <- cbind(SP.dat2,clust.IND) ###
Hy.IND <- hts(Basey[,order(SP.dat2$clust.IND)],nodes=list(nclust.IND,table(clust.IND)))
Ytot.IND <- aggts(Hy.IND)

# 2. Unsupervised approaches:

train.l <- 501 # from 2019 to 2023 training, last year testing

# Returns-based measures:

# EUCL-based clusters (on returns):

Baseret <- apply(Basey, 2, function(x){diff(log(x))})

ASW<-NULL
for (c in 1:4) {
  ASW[c]<-pam(diss(t(Baseret[1:train.l,]), "EUCL"), k=c+1, diss=TRUE)$silinfo$avg.width
}
plot(2:5, ASW, type="b", main="ASW for EUCL-based PAM clustering", ylab="ASW", xlab="Clusters")

nclust.EUCL <- which.max(ASW)+1
clust.EUCL <- pam(diss(t(Baseret[1:train.l,]), "EUCL"), k=nclust.EUCL, diss=TRUE)
Basey.EUCL <- Basey[,order(clust.EUCL$clustering)]
Hy.EUCL <- hts(Basey.EUCL,nodes=list(nclust.EUCL,table(clust.EUCL$clustering)))
Ytot.EUCL <- aggts(Hy.EUCL)

# COR-based clusters (on returns):

ASW<-NULL
for (c in 1:4) {
  ASW[c]<-pam(diss(t(Baseret[1:train.l,]), "COR"), k=c+1, diss=TRUE)$silinfo$avg.width
}

plot(2:5, ASW, type="b", main="ASW for COR-based PAM clustering", ylab="ASW", xlab="Clusters")

nclust.COR <- which.max(ASW)+1
clust.COR<-pam(diss(t(Baseret[1:train.l,]), "COR"), k=nclust.COR, diss=TRUE)
Basey.COR <- Basey[,order(clust.COR$clustering)]
Hy.COR <- hts(Basey.COR,nodes=list(nclust.COR,table(clust.COR$clustering)))
Ytot.COR <- aggts(Hy.COR)

# Prices-based measures:

# ARIMA-based clustering (prices)

cl<-makeCluster(detectCores()-2)
registerDoParallel(cl)

armodels <- foreach(p=1:ncol(Basey), .packages = c("forecast")) %dopar% {
  model0 <- auto.arima(Basey[1:train.l,p], stepwise=F)
  return(model0)
}
names(armodels) <- colnames(Basey)

ModelsAR <- matrix(NA, nrow=ncol(Basey), ncol=1)
for (j in 1:ncol(Basey)) {
  ModelsAR[j,1] <- paste0("ARIMA(",armodels[[j]]$arma[1],",",1,",",armodels[[j]]$arma[2],")")
}
rownames(ModelsAR) <- colnames(Basey)
xtable::xtable(ModelsAR)

k <- 10
PIcoeff <- matrix(NA, k, ncol(Basey))
for (j in 1:ncol(PIcoeff)) {
  model0 <- armodels[[j]]
  PIcoeff[,j] <- ARMAtoAR(ar=model0$model$phi, ma=model0$model$theta, lag.max=k)
}
colnames(PIcoeff)<-colnames(Basey)

plot(cmdscale(dist(t(PIcoeff))))

which(colSums(PIcoeff)==0)# RW

ASW<-NULL
for (c in 1:4) {
  ASW[c]<-pam(diss(t(PIcoeff[-which(colSums(PIcoeff)==0),]), "EUCL"), k=c+1, diss=TRUE)$silinfo$avg.width
}
plot(2:5, ASW, type="b", main="ASW for ARIMA-based PAM clustering", ylab="ASW", xlab="Clusters")

nclust.ARMA <- which.max(ASW)+1
clust.ARMA <- pam(diss(t(PIcoeff[,-which(colSums(PIcoeff)==0)]), "EUCL"), k=nclust.ARMA, diss=TRUE)

Clu3 <- rep(3,length( which(colSums(PIcoeff)==0)))
names(Clu3) <- names(which(colSums(PIcoeff)==0))
CluARMA <- c(clust.ARMA$clustering,Clu3)
clust.ARMA$clustering <- CluARMA[sort(names(CluARMA))]

plot(PIcoeff[,clust.ARMA$clustering==1][,1],ylim=c(min(PIcoeff),max(PIcoeff)),type="l", main="ARMA-based clustering: AR-weights", ylab="Coeff.", xlab="K")
for (i in 2:table(clust.ARMA$clustering)[1]) {
  lines(PIcoeff[,clust.ARMA$clustering==1][,i], col="black")
}
for (i in 1:table(clust.ARMA$clustering)[2]) {
  lines(PIcoeff[,clust.ARMA$clustering==2][,i], col="red")
}
for (i in 1:table(clust.ARMA$clustering)[3]) {
  lines(PIcoeff[,clust.ARMA$clustering==3][,i], col="blue")
}
legend("right",c("Cluster 1", "Cluster 2", "Cluster 3"),col=c("black","red","blue"), lty=1, cex=0.6)

Basey.ARMA <- Basey[,order(clust.ARMA$clustering)]
Hy.ARMA <- hts(Basey.ARMA,nodes=list(3,table(clust.ARMA$clustering)))
Ytot.ARMA <- aggts(Hy.ARMA)

# Combination approaches: ALL

SP.dat2 <- SP.dat[order(SP.dat$Symbol),]
SP.dat2 <- cbind(SP.dat2,clust.IND,clust.EUCL$clustering,
                 clust.COR$clustering, clust.ARMA$clustering)
nclust.all <- nclust.IND+nclust.EUCL+nclust.COR+nclust.ARMA+1

Cmat<-matrix(NA, ncol=ncol(Basey), nrow=nclust.all+1)
colnames(Cmat)<-colnames(Basey)
Cmat[1,]<-rep(1,ncol(Cmat))
for (c in 1:11) {
  Cmat[1+c,]<-ifelse(SP.dat2$clust.IND==c,1,0)
}
for (c in 1:2) {
  Cmat[12+c,]<-ifelse(SP.dat2$`clust.EUCL$clustering`==c,1,0)
}
for (c in 1:2) {
  Cmat[14+c,]<-ifelse(SP.dat2$`clust.COR$clustering`==c,1,0)
}
for (c in 1:3) {
  Cmat[16+c,]<-ifelse(SP.dat2$`clust.ARMA$clustering`==c,1,0)
}
rownames(Cmat)<-1:nrow(Cmat)
rownames(Cmat)[1]<-"tot"
for (c in 1:11) {
  rownames(Cmat)[1+c]<-paste0("c",c,".IND")
}
for (c in 1:2) {
  rownames(Cmat)[12+c]<-paste0("c",c,".EUCL")
}
for (c in 1:2) {
  rownames(Cmat)[14+c]<-paste0("c",c,".COR")
}
for (c in 1:3) {
  rownames(Cmat)[16+c]<-paste0("c",c,".ARMA")
}

Smat <- rbind(Cmat,diag(1,ncol(Basey)))
rownames(Smat)[(nclust.all+2):nrow(Smat)]<-colnames(Basey)
Ytot.all <- matrix(NA, nrow=nrow(Basey), ncol=nrow(Smat))
for (t in 1:nrow(Basey)) {
  Ytot.all[t,] <-  Smat%*%t(Basey[t,])
}
Ytot.tot <- Ytot.all

### B: Forecasting reconciliation

train.l <- 501
out.l <- nrow(Basey)-train.l
h <- c(1,3,6,12)
#h <- c(1,5,22)

Ylists <- list(Ytot.I, Ytot.IND, Ytot.EUCL, Ytot.COR, Ytot.ARMA)
Hlists <- list("NOCLUST"=Hy.I,"IND"=Hy.IND, "EUCL"=Hy.EUCL,"COR"=Hy.COR, "ARMA"=Hy.ARMA)

# Forecasting exercises begins:

cl<-makeCluster(detectCores()-2)
registerDoParallel(cl)

# Set the RDS folder where forecasts will be saved:

setwd(here("RDSforecast_sp500"))

# Cluster-based reconciliation: alternative approaches

for (l in 1:length(h)) {
  
  for (m in 1:5){
    
    Ytot <- Ylists[[m]]
    Hy <- Hlists[[m]] 
    
    #armamodels <- foreach(p=1:ncol(Ytot), .packages = c("forecast")) %dopar% {
    #  model0 <- auto.arima(Ytot[1:train.l,p], seasonal=F, stepwise = F, ic="aic")
    #  return(model0)
    #}
    
    for.h<-foreach(i=1:(out.l-h[l]), .packages = c("forecast")) %dopar% {
      
      yf.h <- matrix(NA, h[l], ncol(Ytot))
      resi <- matrix(NA, train.l, ncol(Ytot)) # train.l+i-1 and Ytot[1:(399+i),j] for expanding window
      
      # ARIMA-based forecasts:
      
      for (j in 1:ncol(Ytot)){
        modelm <- auto.arima(Ytot[i:(500+i),j], stepwise = F)
        yf.h[,j] <- forecast(modelm, h=h[l])$mean#forecast(Arima(Ytot[i:(1203+i),j], model= armamodels[[j]]),h=h[l])$mean 
        resi[,j] <- residuals(modelm) #residuals(Arima(Ytot[i:(1203+i),j], model= armamodels[[j]])) 
      }
      
      return(list(yf.h,resi))
      
    }
    
    export(for.h, paste0("for.h",h[l],names(Hlists)[m],".RDS"))
    
    for.rec<-foreach(i=1:length(for.h), .packages = c("forecast","hts")) %dopar% {
      
      rfor<-MinT(fcasts=for.h[[i]][[1]][h[l],],nodes=get_nodes(Hy), residual=for.h[[i]][[2]], covariance="shr", keep="bottom", algorithms="chol")
      
      return(rfor)
      
    }
    
    export(for.rec, paste0("MinTrec.h",h[l],names(Hlists)[m],".RDS"))
    
    print(paste0(l,"-",m,"completed"))
    
  }
}

#stopCluster(cl)

# ALL reconciliation:

MinTrace <- function(forvec, resimat, Smat, Cmat, keep.bottom=TRUE){
  
  Wmatrix <- hts:::shrink.estim(resimat, hts:::lowerD(resimat))[[1]]
  
  Reconc.f <- Smat%*%solve(t(Smat)%*%solve(Wmatrix)%*%Smat)%*%t(Smat)%*%solve(Wmatrix)%*%t(forvec)
  
  if (keep.bottom == TRUE) {
    
    Reconc.b <- Reconc.f[(nrow(Cmat)+1):nrow(Smat),]
    names(Reconc.b)<-rownames(Smat)[(nrow(Cmat)+1):nrow(Smat)]
    return(Reconc.b)
    
  } else {
    
    return(Reconc.f)
    
  }
} 

for (l in 1:length(h)) {
  
  Ytot0 <- Ytot.all
  Cmat0 <- Cmat
  Smat0 <- Smat
  
  # armamodels <- foreach(p=1:ncol(Ytot0), .packages = c("forecast")) %dopar% {
  #  model0 <- auto.arima(Ytot0[1:train.l,p], seasonal=F, stepwise = F, ic="aic")
  #  return(model0)
  #}
  
  for.h<-foreach(i=1:(out.l-h[l]), .packages = c("forecast")) %dopar% {
    
    yf.h <- matrix(NA, h[l], ncol(Ytot0))
    resi <- matrix(NA, train.l, ncol(Ytot0)) # train.l+i-1 and Ytot0[1:(399+i),j] for expanding window
    
    # ARIMA-based forecasts:
    
    for (j in 1:ncol(Ytot0)) {
      modelm <- auto.arima(Ytot0[i:(500+i),j], stepwise = F)
      yf.h[,j] <- forecast(modelm, h=h[l])$mean#forecast(Arima(Ytot[i:(1203+i),j], model= armamodels[[j]]),h=h[l])$mean 
      resi[,j] <- residuals(modelm) #residuals(Arima(Ytot[i:(1203+i),j], model= armamodels[[j]])) 
    }
    
    return(list(yf.h,resi))
    
  }
  
  export(for.h, paste0("for.h",h[l],"ALL",".RDS"))
  
  for.rec<-foreach(i=1:length(for.h), .packages = c("forecast","hts")) %dopar% {
    
    rfor<-MinTrace(for.h[[i]][[1]], for.h[[i]][[2]], Smat0, Cmat0, keep.bottom=TRUE)
    
    return(rfor)
    
  }
  
  export(for.rec, paste0("MinTrec.h",h[l],"ALL",".RDS"))
  
  print(paste0(l,"completed"))
  
}


# Bottom-up:

for (l in 1:4) {
  
  Base.f <- import(paste0("for.h",h[l],"NOCLUST.RDS"))
  
  Sm <- smatrix(Hy.I)
  G <- cbind(rep(0, ncol(Sm)),diag(1, nrow=ncol(Sm), ncol=ncol(Sm)))
  
  for.rec <- foreach(i=1:length(Base.f), .packages = c("forecast","hts")) %dopar% {
    rfor <- Sm%*%G%*%t(Base.f[[i]][[1]])
    return(rfor)
  }
  
  
  export(for.rec, paste0("rec.h",h[l],"BU",".RDS"))
  
}

# Random walk forecasts:

for (l in 1:length(h)) {
  
  Ytot <- Ytot.I
  
  for.h <- foreach(i = 1:(out.l - h[l]), .packages = c("forecast")) %dopar% {
    
    yf.h <- matrix(NA, h[l], ncol(Ytot))
    resi <- matrix(NA, train.l - 1, ncol(Ytot)) # Use this for RW
    
    # RW-based forecasts:
    for (j in 1:ncol(Ytot)) {
      yf.h[, j] <- rwf(Ytot[i:(500 + i), j], drift = TRUE, h = h[l])$mean
      resi[, j] <- na.omit(resid(rwf(Ytot[i:(500 + i), j], drift = TRUE)))
    }
    
    return(list(yf.h, resi))
  }
  
  saveRDS(for.h, paste0("forRW.h", h[l], ".RDS"))
  
  print(paste0(l, " completed"))
}

stopCluster(cl)

Ytotlist <- list("NOCLUST"=Ytot.I, "IND"=Ytot.IND, "EUCL"=Ytot.EUCL, "COR"=Ytot.COR, 
                 "ARMA"=Ytot.ARMA,"ALL"=Ytot.all)

Hylist <- list("NOCLUST"=Hy.I, "IND"=Hy.IND, "EUCL"=Hy.EUCL, "COR"=Hy.COR, 
               "ARMA"=Hy.ARMA)

export(Ytotlist, "Ytotlist_sp500.RDS")
export(Hylist, "Hylist_sp500.RDS")


########### C. Accuracy evaluation

# Load cluter structure first (part A)

library(rio)
library(here)
library(sandwich)
library(lmtest)

setwd(here("RDSForecast_sp500"))

Ytotlist <- import("Ytotlist_sp500.RDS")

Ytot.I <- Ytotlist[[1]]
Ytot.IND <- Ytotlist[[2]]
Ytot.EUCL <- Ytotlist[[3]]
Ytot.COR <- Ytotlist[[4]]
Ytot.ARMA <- Ytotlist[[5]]
Ytot.all <- Ytotlist[[6]]

train.l <- 501
out.l <- nrow(Ytot.all)-train.l
h <- c(1,3,6,12) # Setting l=1 we get results for h=1 step-ahead; setting h=2 we get h=3 step-ahead, etc.

type <- c("for.h","MinTrec.h")

# Create empty tables

TABtot_index <- matrix(NA, 9, 4) # Tab. 8 (Panel A)
TABtot2_index <- matrix(NA, 9, 4) # Tab. 8 (Panel B)

rownames(TABtot_index) <- c("Base","RW","BU","NO CLUST","IND","EUCL","COR","ARMA","ALL")
rownames(TABtot2_index) <- rownames(TABtot_index)

TABtot_bottom <- matrix(NA, 8, 4) # Tab. 11 (Panel A)
TABtot2_bottom <- matrix(NA, 8, 4) # Tab. 11 (Panel B)

rownames(TABtot_bottom) <- c("Base","RW","NO CLUST","IND","EUCL","COR","ARMA","ALL")
rownames(TABtot2_bottom) <- rownames(TABtot_bottom)

TABpval <- list() # Tab. 9 (each element of the list is a Panel)
TABpval2 <- list() # Tab. 10 (each element of the list is a Panel)
TABpval.bottom <- list() # Tab. 12 (each element of the list is a Panel)
TABpval2.bottom <- list() # Tab. 13 (each element of the list is a Panel)

for (l in 1:4){
  
  print(paste0("Computing",l,"/4..."))
  
  IND.f<-import(paste0(type[1],h[l],"IND.RDS"))
  IND.rec<-import(paste0(type[2],h[l],"IND.RDS"))
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
  
  nclust.IND <- 11 # number of clusters
  nclust.EUCL <- 2 # number of clusters
  nclust.COR <- 2 # number of clusters
  nclust.ARMA <- 3 # number of clusters
  nclust.all <- 18 # number of clusters
  
  ################## DOW JONES INDEX FORECASTING (TOP SERIES)
  
  # 0) BASE WITH RW:
  
  rw.err<-vector("numeric",length(IND.f))
  for (i in 1:length(rw.err)) {
    rw.err[i]<-Ytot.I[(train.l+i+h[l]-1),1]-RW.f[[i]][[1]][h[l],1]
  }
  
  mse.rw<-mean(rw.err^2)
  rmse.rw<-sqrt(mean(rw.err^2))
  mae.rw<-mean(abs(rw.err))
  
  # 1) BASE FORECASTS - DJIA INDEX:
  
  index.err<-vector("numeric",length(IND.f))
  for (i in 1:length(index.err)) {
    index.err[i]<-Ytot.I[(train.l+i+h[l]-1),1]-NOCLUST.f[[i]][[1]][h[l],1]
  }
  
  mse.index<-mean(index.err^2)
  rmse.index<-sqrt(mean(index.err^2))
  mae.index<-mean(abs(index.err))
  
  # 2) BOTTOM UP - DJIA INDEX:
  
  indexBU.err<-vector("numeric",length(ALL.f))
  for (i in 1:length(index.err)) {
    indexBU.err[i]<-Ytot.I[(train.l+i+h[l]-1),1]-sum(BU.rec[[i]][-1,h[l]])
  }
  
  mseBU.index<-mean(indexBU.err^2)
  rmseBU.index<-sqrt(mean(indexBU.err^2))
  maeBU.index<-mean(abs(indexBU.err))
  
  # 3) NO CLUSTERING - DJIA INDEX:
  
  indexNOCLUST.err<-vector("numeric",length(ALL.f))
  for (i in 1:length(index.err)) {
    indexNOCLUST.err[i]<-Ytot.I[(train.l+i+h[l]-1),1]-sum(NOCLUST.rec[[i]][1,])
  }
  
  mseNOCLUST.index<-mean(indexNOCLUST.err^2)
  rmseNOCLUST.index<-sqrt(mean(indexNOCLUST.err^2))
  maeNOCLUST.index<-mean(abs(indexNOCLUST.err))
  
  # 4) IND FORECASTS - DJIA INDEX:
  
  indexIND.err<-vector("numeric",length(IND.f))
  for (i in 1:length(index.err)) {
    indexIND.err[i]<-Ytot.all[(train.l+i+h[l]-1),1]-sum(IND.rec[[i]][1,])
  }
  
  mseIND.index<-mean(indexIND.err^2)
  rmseIND.index<-sqrt(mean(indexIND.err^2))
  maeIND.index<-mean(abs(indexIND.err))
  
  # 6) EUCL FORECASTS - DJIA INDEX:
  
  indexEUCL.err<-vector("numeric",length(EUCL.f))
  for (i in 1:length(index.err)) {
    indexEUCL.err[i]<-Ytot.all[(train.l+i+h[l]-1),1]-sum(EUCL.rec[[i]][1,])
  }
  
  mseEUCL.index<-mean(indexEUCL.err^2)
  rmseEUCL.index<-sqrt(mean(indexEUCL.err^2))
  maeEUCL.index<-mean(abs(indexEUCL.err))
  
  # 7) COR FORECASTS - DJIA INDEX:
  
  indexCOR.err<-vector("numeric",length(COR.f))
  for (i in 1:length(index.err)) {
    indexCOR.err[i]<-Ytot.all[(train.l+i+h[l]-1),1]-sum(COR.rec[[i]][1,])
  }
  
  mseCOR.index<-mean(indexCOR.err^2)
  rmseCOR.index<-sqrt(mean(indexCOR.err^2))
  maeCOR.index<-mean(abs(indexCOR.err))
  
  # 8) ARMA FORECASTS - DJIA INDEX:
  
  indexARMA.err<-vector("numeric",length(ARMA.f))
  for (i in 1:length(index.err)) {
    indexARMA.err[i]<-Ytot.all[(train.l+i+h[l]-1),1]-sum(ARMA.rec[[i]][1,])
  }
  
  mseARMA.index<-mean(indexARMA.err^2)
  rmseARMA.index<-sqrt(mean(indexARMA.err^2))
  maeARMA.index<-mean(abs(indexARMA.err))
  
  # 9) ALL FORECASTS - DJIA INDEX:
  
  indexALL.err<-vector("numeric",length(ALL.f))
  for (i in 1:length(index.err)) {
    indexALL.err[i]<-ifelse(h[l]==1, Ytot.all[(train.l+i+h[l]-1),1]-sum(ALL.rec[[i]]), Ytot.all[(train.l+i+h[l]-1),1]-sum(ALL.rec[[i]][,h[l]]))
  }
  
  mseALL.index<-mean(indexALL.err^2)
  rmseALL.index<-sqrt(mean(indexALL.err^2))
  maeALL.index<-mean(abs(indexALL.err))
  
  # Comparisons:
  
  ERRORS.index<-cbind(index.err,rw.err,indexBU.err,indexNOCLUST.err,indexIND.err,indexEUCL.err,indexCOR.err,indexARMA.err,indexALL.err)
  
  MAE <- abs(ERRORS.index)
  MSE <- ERRORS.index^2
  
  # Compute MAE metric:
  
  TABtot_index[,l] <- colMeans(MAE)
  
  
  # Compute RMSE metric: 
  
  TABtot2_index[,l] <- sqrt(colMeans(MSE))
  
  # Take function from the empirical density of Y:
  
  Y0t <- Ytot.all[(train.l+1+h[l]):nrow(Ytot.all),1]
  d <- density(Y0t)
  density_fun <- approxfun(d$x, d$y, method="constant")
  w <- density_fun(Y0t)
  w <- (w-min(w))/(max(w)-min(w))
  
  # Inference with absolute errors:
  
  TABpval[[l]] <- matrix(NA, nrow=5, ncol=4)
  
  for (i in 1:nrow(TABpval[[l]])) {
    derr <- MAE[,1]-MAE[,4+i]
    derr <- w*derr
    lm0 <- lm(derr~1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
    TABpval[[l]][i, 1] <- pt(abs(derr.test[1]/derr.test[2]), df=length(derr)-1, lower.tail = FALSE)
    
    derr <- MAE[,2]-MAE[,4+i]
    derr <- w*derr
    lm0 <- lm(derr~1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
    TABpval[[l]][i, 2] <- pt(abs(derr.test[1]/derr.test[2]), df=length(derr)-1, lower.tail = FALSE)
    
    derr <- MAE[,3]-MAE[,4+i]
    derr <- w*derr
    lm0 <- lm(derr~1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
    TABpval[[l]][i, 3] <- pt(abs(derr.test[1]/derr.test[2]), df=length(derr)-1, lower.tail = FALSE)
    
    derr <- MAE[,4]-MAE[,4+i]
    derr <- w*derr
    lm0 <- lm(derr~1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
    TABpval[[l]][i, 4] <- pt(abs(derr.test[1]/derr.test[2]), df=length(derr)-1, lower.tail = FALSE)
  }
  
  colnames(TABpval[[l]]) <- c("Base","RW","BU","NO CLUST")
  rownames(TABpval[[l]]) <- c("IND","EUCL","COR","ARMA","ALL")
  
  # Inference with squared errors:
  
  TABpval2[[l]] <-matrix(NA, nrow=5, ncol=4)
  
  for (i in 1:nrow(TABpval2[[l]])) {
    derr <- MSE[,1]-MSE[,4+i]
    derr <- w*derr
    lm0 <- lm(derr~1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
    TABpval2[[l]][i, 1] <- pt(abs(derr.test[1]/derr.test[2]), df = length(derr)-1, lower.tail = FALSE)
    
    derr <- MSE[,2]-MSE[,4+i]
    derr <- w*derr
    lm0 <- lm(derr~1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
    TABpval2[[l]][i, 2] <- pt(abs(derr.test[1]/derr.test[2]), df = length(derr)-1, lower.tail = FALSE)
    
    derr <- MSE[,3]-MSE[,4+i]
    derr <- w*derr
    lm0 <- lm(derr~1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
    TABpval2[[l]][i, 3] <- pt(abs(derr.test[1]/derr.test[2]), df = length(derr)-1, lower.tail = FALSE)
    
    derr <- MSE[,4]-MSE[,4+i]
    derr <- w*derr
    lm0 <- lm(derr~1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
    TABpval2[[l]][i, 4] <- pt(abs(derr.test[1]/derr.test[2]), df = length(derr)-1, lower.tail = FALSE)
  }
  
  colnames(TABpval2[[l]]) <- c("Base","RW","BU","NO CLUST")
  rownames(TABpval2[[l]]) <- c("IND","EUCL","COR","ARMA","ALL")
  
  ################## COSTITUENTS FORECASTING (BOTTOM SERIES)
  
  # 1) BASE FORECASTS:
  
  # Are the same of bottom-up:
  
  fbase.err<-matrix(NA, nrow=length(IND.f), ncol=497)
  fbase.perr<-matrix(NA, nrow=length(IND.f), ncol=497)
  for (i in 1:nrow(fbase.err)) {
    fbase.err[i,]<-Ytot.I[(train.l+i+h[l]-1),2:ncol(Ytot.I)]-BU.rec[[i]][-1,h[l]] # BU.rec -1 to remove the top series
  }
  for (i in 1:nrow(fbase.err)) {
    fbase.perr[i,]<-100*(Ytot.I[(train.l+i+h[l]-1),2:ncol(Ytot.I)]-BU.rec[[i]][-1,h[l]])/Ytot.I[(train.l+i+h[l]-1),2:ncol(Ytot.I)]
  }
  mspe.base<-apply(fbase.err, 2, function(x){sqrt(mean(x^2))}) #log(x^2+1)
  mape.base<-apply(fbase.err, 2, function(x){mean(abs(x))}) # log(abs(x+1))
  mean(mape.base)
  mean(mspe.base)
  
  # 2) RW:
  
  RW.err<-matrix(NA, nrow=length(IND.f), ncol=497)
  RW.perr<-matrix(NA, nrow=length(IND.f), ncol=497)
  for (i in 1:nrow(RW.err)) {
    RW.err[i,]<-Ytot.I[(train.l+i+h[l]-1),-1]-RW.f[[i]][[1]][h[l],-1]
  }
  for (i in 1:nrow(RW.err)) {
    RW.perr[i,]<-100*(Ytot.I[(train.l+i+h[l]-1),-1]-RW.f[[i]][[1]][h[l],-1])/Ytot.I[(train.l+i+h[l]-1),-1]
  }
  mspe.RW<-apply(RW.err, 2, function(x){sqrt(mean(x^2))})
  mape.RW<-apply(RW.err, 2, function(x){mean(abs(x))})
  mean(mape.RW)
  mean(mspe.RW)
  
  # 3) No clustering
  
  recNOCLUST.err<-matrix(NA, nrow=length(ALL.f), ncol=497)
  recNOCLUST.perr<-matrix(NA, nrow=length(ALL.f), ncol=497)
  for (i in 1:nrow(fbase.err)) {
    recNOCLUST.err[i,]<-Ytot.I[(train.l+i+h[l]-1),2:ncol(Ytot.I)]-NOCLUST.rec[[i]][1,] # ALL.rec[[i]] if h=1, ALL.rec[[i]][,h[l]] for h>1.
  }
  for (i in 1:nrow(fbase.err)) {
    recNOCLUST.perr[i,]<-100*(Ytot.I[(train.l+i+h[l]-1),2:ncol(Ytot.I)]-NOCLUST.rec[[i]][1,])/Ytot.I[(train.l+i+h[l]-1),2:ncol(Ytot.I)]
  }
  mspeNOCLUST.base<-apply(recNOCLUST.err, 2, function(x){sqrt(mean(x^2))})
  mapeNOCLUST.base<-apply(recNOCLUST.err, 2, function(x){mean(abs(x))})
  mean(mapeNOCLUST.base)
  mean(mspeNOCLUST.base)
  
  # 4) IND-based reconcilied forecasts:
  
  recIND.err<-matrix(NA, nrow=length(IND.f), ncol=497)
  recIND.perr<-matrix(NA, nrow=length(IND.f), ncol=497)
  for (i in 1:nrow(fbase.err)) {
    recIND.err[i,]<-Ytot.IND[(train.l+i+h[l]-1),(nclust.IND+2):ncol(Ytot.IND)]-IND.rec[[i]][1,]
  }
  for (i in 1:nrow(fbase.err)) {
    recIND.perr[i,]<-100*(Ytot.IND[(train.l+i+h[l]-1),(nclust.IND+2):ncol(Ytot.IND)]-IND.rec[[i]][1,])/Ytot.IND[(train.l+i+h[l]-1),(nclust.IND+2):ncol(Ytot.IND)]
}
mspeIND.base<-apply(recIND.err, 2, function(x){sqrt(mean(x^2))})
mapeIND.base<-apply(recIND.err, 2, function(x){mean(abs(x))})
mean(mapeIND.base)
mean(mspeIND.base)

# 5) EUCL-based reconcilied forecasts:

recEUCL.err<-matrix(NA, nrow=length(EUCL.f), ncol=497)
recEUCL.perr<-matrix(NA, nrow=length(EUCL.f), ncol=497)
for (i in 1:nrow(fbase.err)) {
  recEUCL.err[i,]<-Ytot.EUCL[(train.l+i+h[l]-1),(nclust.EUCL+2):ncol(Ytot.EUCL)]-EUCL.rec[[i]][1,]
}
for (i in 1:nrow(fbase.err)) {
  recEUCL.perr[i,]<-100*(Ytot.EUCL[(train.l+i+h[l]-1),(nclust.EUCL+2):ncol(Ytot.EUCL)]-EUCL.rec[[i]][1,])/Ytot.EUCL[(train.l+i+h[l]-1),(nclust.EUCL+2):ncol(Ytot.EUCL)]
}
mspeEUCL.base<-apply(recEUCL.err, 2, function(x){sqrt(mean(x^2))})
mapeEUCL.base<-apply(recEUCL.err, 2, function(x){mean(abs(x))})
mean(mapeEUCL.base)
mean(mspeEUCL.base)

# 7) COR-based reconcilied forecasts:

recCOR.err<-matrix(NA, nrow=length(COR.f), ncol=497)
recCOR.perr<-matrix(NA, nrow=length(COR.f), ncol=497)
for (i in 1:nrow(fbase.err)) {
  recCOR.err[i,]<-Ytot.COR[(train.l+i+h[l]-1),(nclust.COR+2):ncol(Ytot.COR)]-COR.rec[[i]][1,]
}
for (i in 1:nrow(fbase.err)) {
  recCOR.perr[i,]<-100*(Ytot.COR[(train.l+i+h[l]-1),(nclust.COR+2):ncol(Ytot.COR)]-COR.rec[[i]][1,])/Ytot.COR[(train.l+i+h[l]-1),(nclust.COR+2):ncol(Ytot.COR)]
}
mspeCOR.base<-apply(recCOR.err, 2, function(x){sqrt(mean(x^2))})
mapeCOR.base<-apply(recCOR.err, 2, function(x){mean(abs(x))})
mean(mapeCOR.base)
mean(mspeCOR.base)

# 8) ARMA-based reconcilied forecasts:

recARMA.err<-matrix(NA, nrow=length(ARMA.f), ncol=497)
recARMA.perr<-matrix(NA, nrow=length(ARMA.f), ncol=497)
for (i in 1:nrow(fbase.err)) {
  recARMA.err[i,]<-Ytot.ARMA[(train.l+i+h[l]-1),(nclust.ARMA+2):ncol(Ytot.ARMA)]-ARMA.rec[[i]][1,]
}
for (i in 1:nrow(fbase.err)) {
  recARMA.perr[i,]<-100*(Ytot.ARMA[(train.l+i+h[l]-1),(nclust.ARMA+2):ncol(Ytot.ARMA)]-ARMA.rec[[i]][1,])/Ytot.ARMA[(train.l+i+h[l]-1),(nclust.ARMA+2):ncol(Ytot.ARMA)]
}
mspeARMA.base<-apply(recARMA.err, 2, function(x){sqrt(mean(x^2))})
mapeARMA.base<-apply(recARMA.err, 2, function(x){mean(abs(x))})
mean(mapeARMA.base)
mean(mspeARMA.base)

# 9) ALL-based reconcilied forecasts:

recALL.err<-matrix(NA, nrow=length(ALL.f), ncol=497)
recALL.perr<-matrix(NA, nrow=length(ALL.f), ncol=497)

if (h[l]==1){
for (i in 1:nrow(fbase.err)) {
  recALL.err[i,] <- Ytot.all[(train.l+i+h[l]-1),(nclust.all+2):ncol(Ytot.all)]-ALL.rec[[i]]
}
for (i in 1:nrow(fbase.err)) {
  recALL.perr[i,] <- 100*(Ytot.all[(train.l+i+h[l]-1),(nclust.all+2):ncol(Ytot.all)]-ALL.rec[[i]])/Ytot.all[(train.l+i+h[l]-1),(nclust.all+2):ncol(Ytot.all)]
} } else {
  for (i in 1:nrow(fbase.err)) {
    recALL.err[i,]<- Ytot.all[(train.l+i+h[l]-1),(nclust.all+2):ncol(Ytot.all)]-ALL.rec[[i]][,h[l]] # ALL.rec[[i]] if h=1, ALL.rec[[i]][,h[l]] for h>1.
  }
  for (i in 1:nrow(fbase.err)) {
    recALL.perr[i,]<- 100*(Ytot.all[(train.l+i+h[l]-1),(nclust.all+2):ncol(Ytot.all)]-ALL.rec[[i]][,h[l]])/Ytot.all[(train.l+i+h[l]-1),(nclust.all+2):ncol(Ytot.all)]
  }  
}
mspeALL.base<-apply(recALL.err, 2, function(x){sqrt(mean(x^2))})
mapeALL.base<-apply(recALL.err, 2, function(x){mean(abs(x))})
mean(mapeALL.base)
mean(mspeALL.base)

# Accuracy measures comparison - average

# MAE:


MAE <- cbind(rowMeans(abs(fbase.perr)),rowMeans(abs(RW.perr)),rowMeans(abs(recNOCLUST.perr)),
             rowMeans(abs(recIND.perr)), rowMeans(abs(recEUCL.perr)),
             rowMeans(abs(recCOR.perr)), rowMeans(abs(recARMA.perr)), rowMeans(abs(recALL.perr)))

MAE.BASE<-rbind(mean(mape.base), mean(mape.RW), mean(mapeNOCLUST.base),
                mean(mapeIND.base), 
                mean(mapeEUCL.base), mean(mapeCOR.base),
                mean(mapeARMA.base), mean(mapeALL.base))


names(MAE.BASE)<-c("Base","RW","NO CLUST","IND","EUCL","COR","ARMA","ALL")

TABtot_bottom[,l] <- MAE.BASE

# MSE:

MSE <- cbind(rowMeans(fbase.perr^2),rowMeans(RW.perr^2),rowMeans(recNOCLUST.perr^2),
             rowMeans(recIND.perr^2), rowMeans(recEUCL.perr^2),
             rowMeans(recCOR.perr^2), rowMeans(recARMA.perr^2), rowMeans(recALL.perr^2))

MSE.BASE<-rbind(mean(mspe.base), mean(mspe.RW), mean(mspeNOCLUST.base),
                mean(mspeIND.base),  
                mean(mspeEUCL.base), mean(mspeCOR.base),
                mean(mspeARMA.base), mean(mspeALL.base))

names(MSE.BASE)<-c("Base","RW","NO CLUST","IND","EUCL","COR","ARMA","ALL")

TABtot2_bottom[,l] <- MSE.BASE

# Take function from the empirical density of Y:

Y0t <- Ytot.all[(train.l+1+h[l]):nrow(Ytot.all),1]
d <- density(Y0t)
density_fun <- approxfun(d$x, d$y, method="constant")
w <- density_fun(Y0t)
w <- (w-min(w))/(max(w))

# Inference with MAE loss:

TABpval.bottom[[l]] <-matrix(NA, nrow=5, ncol=3)

for (i in 1:nrow(TABpval.bottom[[l]])) {
  derr <- MAE[,1]-MAE[,3+i]
  #derr <- w*derr
  lm0 <- lm(derr~1)
  derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
  TABpval.bottom[[l]][i, 1] <- pt(abs(derr.test[1]/derr.test[2]), df = length(derr)-1, lower.tail = FALSE)
  
  derr <- MAE[,2]-MAE[,3+i]
  #derr <- w*derr
  lm0 <- lm(derr~1)
  derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
  TABpval.bottom[[l]][i, 2] <- pt(abs(derr.test[1]/derr.test[2]), df = length(derr)-1, lower.tail = FALSE)
  
  derr <- MAE[,3]-MAE[,3+i]
  #derr <- w*derr
  lm0 <- lm(derr~1)
  derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
  TABpval.bottom[[l]][i, 3] <- pt(abs(derr.test[1]/derr.test[2]), df = length(derr)-1, lower.tail = FALSE)
  
}

colnames(TABpval.bottom[[l]]) <- c("Base","RW","NO CLUST")
rownames(TABpval.bottom[[l]]) <- c("IND","EUCL","COR","ARMA","ALL")

# Inference with RMSE loss:

TABpval2.bottom[[l]] <-matrix(NA, nrow=5, ncol=3)

for (i in 1:nrow(TABpval2.bottom[[l]])) {
  derr <- MSE[,1]-MSE[,3+i]
  #derr <- w*derr
  lm0 <- lm(derr~1)
  derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
  TABpval2.bottom[[l]][i, 1] <- pt(abs(derr.test[1]/derr.test[2]), df = length(derr)-1, lower.tail = FALSE)
  
  derr <- MSE[,2]-MSE[,3+i]
  #derr <- w*derr
  lm0 <- lm(derr~1)
  derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
  TABpval2.bottom[[l]][i, 2] <- pt(abs(derr.test[1]/derr.test[2]), df = length(derr)-1, lower.tail = FALSE)
  
  derr <- MSE[,3]-MSE[,3+i]
  #derr <- w*derr
  lm0 <- lm(derr~1)
  derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag=100))
  TABpval2.bottom[[l]][i, 3] <- pt(abs(derr.test[1]/derr.test[2]), df = length(derr)-1, lower.tail = FALSE)
  
}

colnames(TABpval2.bottom[[l]]) <- c("Base","RW","NO CLUST")
rownames(TABpval2.bottom[[l]]) <- c("IND","EUCL","COR","ARMA","ALL")

}

TABtot_index
TABtot2_index

xtable::xtable(TABtot_index, digits=2)
xtable::xtable(TABtot2_index, digits=2)

TABtot_bottom
TABtot2_bottom

xtable::xtable(TABtot_bottom, digits=2)
xtable::xtable(TABtot2_bottom, digits=2)

TABpval
TABpval2

TABpval.bottom
TABpval2.bottom

########### Timing

SR.test <- function(Price, predPrice, pred2, c=0.005) {
  
  # Compute predicted (and actual) returns for direction:
  
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
  
  # Compute net returns given cost c
  
  cost <- NULL
  for (i in 1:length(Price)) {
    cost[i] <- (1+c)*Price[i]
  }
  #cost <- na.omit(cost)
  
  ret1 <- NULL
  for (i in 1:length(Price)) {
    ret1[i] <- (Price[i]-cost[i-1])/cost[i-1]
  }
  ret1 <- na.omit(ret1)
  
  ret2 <- NULL
  for (i in 1:length(Price)) {
    ret2[i] <- (Price[i]-cost[i-1])/cost[i-1]
  }
  ret2 <- na.omit(ret2)
  
  r_t = sign1 * ret1
  r_t2 = sign2 * ret2
  
  stat = ((mean(r_t)/sd(r_t)) - (mean(r_t2)/sd(r_t2))) * 100
  Shr = sd(r_t2)*mean(r_t) - sd(r_t)*mean(r_t2)
  Theta = (var(r_t)*var(r_t2)/length(ret)) * ( 2*cov(r_t,r_t2) + 0.5 * ( (mean(r_t)/sd(r_t))^2 + (mean(r_t2)/sd(r_t2))^2 - ((mean(r_t)/sd(r_t))*(mean(r_t2)/sd(r_t2)) * (1+cov(r_t, r_t2)^2)) ) )
  
  p.val = 1-pnorm(Shr/sqrt(Theta))
  
  return(list=c("Statistic"=stat, p.val=p.val))
  
}

Ytot.I <- Ytotlist[[1]]
Ytot.IND <- Ytotlist[[2]]
Ytot.EUCL <- Ytotlist[[3]]
Ytot.COR <- Ytotlist[[4]]
Ytot.ARMA <- Ytotlist[[5]]
Ytot.all <- Ytotlist[[6]]

train.l <- 501
out.l <- nrow(Ytot.all)-train.l
h <- c(1,3,6,12) 
type <- c("for.h","MinTrec.h")

InvRes <- list()
InvPval <- list()

for (l in 1:4){
  
  print(paste0("Computing",l,"/4..."))
  
  IND.f<-import(paste0(type[1],h[l],"IND.RDS"))
  IND.rec<-import(paste0(type[2],h[l],"IND.RDS"))
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
    yhat.IND[i] <- sum(IND.rec[[i]][1,])
  }
  
  yhat.EUCL <- NULL
  for (i in 1:length(IND.f)) {
    yhat.EUCL[i] <- sum(EUCL.rec[[i]][1,])
  }
  
  yhat.COR <- NULL
  for (i in 1:length(IND.f)) {
    yhat.COR[i] <- sum(COR.rec[[i]][1,])
  }
  
  yhat.ARMA <- NULL
  for (i in 1:length(IND.f)) {
    yhat.ARMA[i] <- sum(ARMA.rec[[i]][1,])
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
    yhat.NOCLUST[i] <- sum(NOCLUST.rec[[i]][1,])
  }
  
  yhat.Base <- NULL
  for (i in 1:length(IND.f)) {
    yhat.Base[i] <- NOCLUST.f[[i]][[1]][h[l],1]
  }
  
  
  # Proportion test: ARMA vs reconciliation h=1
  
  InvMat1 <- matrix(NA,5,3)
  InvMat2 <- matrix(NA,5,3)
  
  InvMat1[1,1] <- SR.test(yval, yhat.IND, yhat.Base)[1]
  InvMat2[1,1] <- SR.test(yval, yhat.IND, yhat.Base)[2]
  InvMat1[2,1] <- SR.test(yval, yhat.EUCL, yhat.Base)[1]
  InvMat2[2,1] <- SR.test(yval, yhat.EUCL, yhat.Base)[2]
  InvMat1[3,1] <- SR.test(yval, yhat.COR, yhat.Base)[1]
  InvMat2[3,1] <- SR.test(yval, yhat.COR, yhat.Base)[2]
  InvMat1[4,1] <- SR.test(yval, yhat.ARMA, yhat.Base)[1]
  InvMat2[4,1] <- SR.test(yval, yhat.ARMA, yhat.Base)[2]
  InvMat1[5,1] <- SR.test(yval, yhat.ALL, yhat.Base)[1]
  InvMat2[5,1] <- SR.test(yval, yhat.ALL, yhat.Base)[2]
  
  InvMat1[1,2] <- SR.test(yval, yhat.IND, yhat.BU)[1]
  InvMat2[1,2] <- SR.test(yval, yhat.IND, yhat.BU)[2]
  InvMat1[2,2] <- SR.test(yval, yhat.EUCL, yhat.BU)[1]
  InvMat2[2,2] <- SR.test(yval, yhat.EUCL, yhat.BU)[2]
  InvMat1[3,2] <- SR.test(yval, yhat.COR, yhat.BU)[1]
  InvMat2[3,2] <- SR.test(yval, yhat.COR, yhat.BU)[2]
  InvMat1[4,2] <- SR.test(yval, yhat.ARMA, yhat.BU)[1]
  InvMat2[4,2] <- SR.test(yval, yhat.ARMA, yhat.BU)[2]
  InvMat1[5,2] <- SR.test(yval, yhat.ALL, yhat.BU)[1]
  InvMat2[5,2] <- SR.test(yval, yhat.ALL, yhat.BU)[2]
  
  InvMat1[1,3] <- SR.test(yval, yhat.IND, yhat.NOCLUST)[1]
  InvMat2[1,3] <- SR.test(yval, yhat.IND, yhat.NOCLUST)[2]
  InvMat1[2,3] <- SR.test(yval, yhat.EUCL, yhat.NOCLUST)[1]
  InvMat2[2,3] <- SR.test(yval, yhat.EUCL, yhat.NOCLUST)[2]
  InvMat1[3,3] <- SR.test(yval, yhat.COR, yhat.NOCLUST)[1]
  InvMat2[3,3] <- SR.test(yval, yhat.COR, yhat.NOCLUST)[2]
  InvMat1[4,3] <- SR.test(yval, yhat.ARMA, yhat.NOCLUST)[1]
  InvMat2[4,3] <- SR.test(yval, yhat.ARMA, yhat.NOCLUST)[2]
  InvMat1[5,3] <- SR.test(yval, yhat.ALL, yhat.NOCLUST)[1]
  InvMat2[5,3] <- SR.test(yval, yhat.ALL, yhat.NOCLUST)[2]
  
  InvRes[[l]] <- InvMat1
  InvPval[[l]] <- InvMat2
  
}

print(xtable::xtable(InvRes[[1]], digits=2), include.rownames = F)
print(xtable::xtable(InvRes[[2]], digits=2), include.rownames = F)
round(InvPval[[1]],2)
