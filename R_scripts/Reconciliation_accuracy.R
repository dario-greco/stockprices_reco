#
# Code for reproducing results in Section 5
#
###########################################################

library(rio)
library(here)
library(sandwich)
library(lmtest)
library(tidyverse)

# Default colours
options(
  ggplot2.discrete.colour = c("#D55E00", "#0072B2","#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442"),
  ggplot2.discrete.fill = c("#D55E00", "#0072B2","#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442")
)

Ytotlist <- import(here("RDSForecast/Ytotlist.RDS"))

Ytot.I <- Ytotlist[[1]]
Ytot.IND <- Ytotlist[[2]]
Ytot.MRKT <- Ytotlist[[3]]
Ytot.EUCL <- Ytotlist[[4]]
Ytot.COR <- Ytotlist[[5]]
Ytot.ARMA <- Ytotlist[[6]]
Ytot.all <- Ytotlist[[7]]

train.l <- 400
out.l <- nrow(Ytot.all) - train.l
h <- c(1, 3, 6, 12) # Setting l=1 we get results for h=1 step-ahead; setting h=2 we get h=3 step-ahead, etc.

# Plotting out-of-sample series:

pdf(here("Figures/oosP.pdf"), width = 8, height = 5)
plot(Ytot.I[(train.l + 1):nrow(Ytot.I), 1], type = "l", main = "DJIA price: testing set", ylab = "Price")
crop::dev.off.crop()

type <- c("for.h", "MinTrec.h")

# Create empty tables

TABtot_index <- matrix(NA, 10, 4) # Tab. 8 (Panel A)
TABtot2_index <- matrix(NA, 10, 4) # Tab. 8 (Panel B)

rownames(TABtot_index) <- c("Base", "RW", "BU", "NO CLUST", "IND", "MRKT", "EUCL", "COR", "ARMA", "ALL")
rownames(TABtot2_index) <- rownames(TABtot_index)

TABtot_bottom <- matrix(NA, 9, 4) # Tab. 11 (Panel A)
TABtot2_bottom <- matrix(NA, 9, 4) # Tab. 11 (Panel B)

rownames(TABtot_bottom) <- c("Base", "RW", "NO CLUST", "IND", "MRKT", "EUCL", "COR", "ARMA", "ALL")
rownames(TABtot2_bottom) <- rownames(TABtot_bottom)

TABpval <- list() # Tab. 9 (each element of the list is a Panel)
TABpval2 <- list() # Tab. 10 (each element of the list is a Panel)
TABpval.bottom <- list() # Tab. 12 (each element of the list is a Panel)
TABpval2.bottom <- list() # Tab. 13 (each element of the list is a Panel)

for (l in 1:4) {
  print(paste0("Computing", l, "/4..."))

  IND.f <- import(here(paste0("RDSForecast/", type[1], h[l], "IND.RDS")))
  IND.rec <- import(here(paste0("RDSForecast/", type[2], h[l], "IND.RDS")))
  MRKT.f <- import(here(paste0("RDSForecast/", type[1], h[l], "MRKT.RDS")))
  MRKT.rec <- import(here(paste0("RDSForecast/", type[2], h[l], "MRKT.RDS")))
  EUCL.f <- import(here(paste0("RDSForecast/", type[1], h[l], "EUCL.RDS")))
  EUCL.rec <- import(here(paste0("RDSForecast/", type[2], h[l], "EUCL.RDS")))
  COR.f <- import(here(paste0("RDSForecast/", type[1], h[l], "COR.RDS")))
  COR.rec <- import(here(paste0("RDSForecast/", type[2], h[l], "COR.RDS")))
  ARMA.f <- import(here(paste0("RDSForecast/", type[1], h[l], "ARMA.RDS")))
  ARMA.rec <- import(here(paste0("RDSForecast/", type[2], h[l], "ARMA.RDS")))
  ALL.f <- import(here(paste0("RDSForecast/", type[1], h[l], "ALL.RDS")))
  ALL.rec <- import(here(paste0("RDSForecast/", type[2], h[l], "ALL.RDS")))
  NOCLUST.f <- import(here(paste0("RDSForecast/", type[1], h[l], "NOCLUST.RDS")))
  NOCLUST.rec <- import(here(paste0("RDSForecast/", type[2], h[l], "NOCLUST.RDS")))
  BU.rec <- import(here(paste0("RDSForecast/rec.h", h[l], "BU", ".RDS")))
  RW.f <- import(here(paste0("RDSForecast/forRW.h", h[l], ".RDS")))

  nclust.IND <- 20 # number of clusters
  nclust.MRKT <- 2 # number of clusters
  nclust.EUCL <- 3 # number of clusters
  nclust.COR <- 3 # number of clusters
  nclust.ARMA <- 5 # number of clusters
  nclust.all <- 33 # number of clusters

  ################## DOW JONES INDEX FORECASTING (TOP SERIES)

  # 0) BASE WITH RW:

  rw.err <- vector("numeric", length(IND.f))
  for (i in 1:length(rw.err)) {
    rw.err[i] <- Ytot.I[(train.l + i + h[l] - 1), 1] - RW.f[[i]][[1]][h[l], 1]
  }

  mse.rw <- mean(rw.err^2)
  rmse.rw <- sqrt(mean(rw.err^2))
  mae.rw <- mean(abs(rw.err))

  # 1) BASE FORECASTS - DJIA INDEX:

  index.err <- vector("numeric", length(IND.f))
  for (i in 1:length(index.err)) {
    index.err[i] <- Ytot.I[(train.l + i + h[l] - 1), 1] - NOCLUST.f[[i]][[1]][h[l], 1]
  }

  mse.index <- mean(index.err^2)
  rmse.index <- sqrt(mean(index.err^2))
  mae.index <- mean(abs(index.err))

  # 2) BOTTOM UP - DJIA INDEX:

  indexBU.err <- vector("numeric", length(ALL.f))
  for (i in 1:length(index.err)) {
    indexBU.err[i] <- Ytot.I[(train.l + i + h[l] - 1), 1] - sum(BU.rec[[i]][-1, h[l]])
  }

  mseBU.index <- mean(indexBU.err^2)
  rmseBU.index <- sqrt(mean(indexBU.err^2))
  maeBU.index <- mean(abs(indexBU.err))

  # 3) NO CLUSTERING - DJIA INDEX:

  indexNOCLUST.err <- vector("numeric", length(ALL.f))
  for (i in 1:length(index.err)) {
    indexNOCLUST.err[i] <- Ytot.I[(train.l + i + h[l] - 1), 1] - sum(NOCLUST.rec[[i]][1, ]) # 1 for the clustered approaches (except ALL) means the h[l] step ahead forecast. Ex: instead of all forecasts for 1,2,3 we just take 3
  }

  mseNOCLUST.index <- mean(indexNOCLUST.err^2)
  rmseNOCLUST.index <- sqrt(mean(indexNOCLUST.err^2))
  maeNOCLUST.index <- mean(abs(indexNOCLUST.err))

  # 4) IND FORECASTS - DJIA INDEX:

  indexIND.err <- vector("numeric", length(IND.f))
  for (i in 1:length(index.err)) {
    indexIND.err[i] <- Ytot.all[(train.l + i + h[l] - 1), 1] - sum(IND.rec[[i]][1, ])
  }

  mseIND.index <- mean(indexIND.err^2)
  rmseIND.index <- sqrt(mean(indexIND.err^2))
  maeIND.index <- mean(abs(indexIND.err))

  # 5) MRKT FORECASTS - DJIA INDEX:

  indexMRKT.err <- vector("numeric", length(MRKT.f))
  for (i in 1:length(index.err)) {
    indexMRKT.err[i] <- Ytot.all[(train.l + i + h[l] - 1), 1] - sum(MRKT.rec[[i]][1, ])
  }

  mseMRKT.index <- mean(indexMRKT.err^2)
  rmseMRKT.index <- sqrt(mean(indexMRKT.err^2))
  maeMRKT.index <- mean(abs(indexMRKT.err))

  # 6) EUCL FORECASTS - DJIA INDEX:

  indexEUCL.err <- vector("numeric", length(EUCL.f))
  for (i in 1:length(index.err)) {
    indexEUCL.err[i] <- Ytot.all[(train.l + i + h[l] - 1), 1] - sum(EUCL.rec[[i]][1, ])
  }

  mseEUCL.index <- mean(indexEUCL.err^2)
  rmseEUCL.index <- sqrt(mean(indexEUCL.err^2))
  maeEUCL.index <- mean(abs(indexEUCL.err))

  # 7) COR FORECASTS - DJIA INDEX:

  indexCOR.err <- vector("numeric", length(COR.f))
  for (i in 1:length(index.err)) {
    indexCOR.err[i] <- Ytot.all[(train.l + i + h[l] - 1), 1] - sum(COR.rec[[i]][1, ])
  }

  mseCOR.index <- mean(indexCOR.err^2)
  rmseCOR.index <- sqrt(mean(indexCOR.err^2))
  maeCOR.index <- mean(abs(indexCOR.err))

  # 8) ARMA FORECASTS - DJIA INDEX:

  indexARMA.err <- vector("numeric", length(ARMA.f))
  for (i in 1:length(index.err)) {
    indexARMA.err[i] <- Ytot.all[(train.l + i + h[l] - 1), 1] - sum(ARMA.rec[[i]][1, ])
  }

  mseARMA.index <- mean(indexARMA.err^2)
  rmseARMA.index <- sqrt(mean(indexARMA.err^2))
  maeARMA.index <- mean(abs(indexARMA.err))

  # 9) ALL FORECASTS - DJIA INDEX:

  indexALL.err <- vector("numeric", length(ALL.f))
  for (i in 1:length(index.err)) {
    indexALL.err[i] <- ifelse(h[l] == 1, Ytot.all[(train.l + i + h[l] - 1), 1] - sum(ALL.rec[[i]]), Ytot.all[(train.l + i + h[l] - 1), 1] - sum(ALL.rec[[i]][, h[l]]))
  }

  mseALL.index <- mean(indexALL.err^2)
  rmseALL.index <- sqrt(mean(indexALL.err^2))
  maeALL.index <- mean(abs(indexALL.err))

  # Comparisons:

  ERRORS.index <- cbind(index.err, rw.err, indexBU.err, indexNOCLUST.err, indexIND.err, indexMRKT.err, indexEUCL.err, indexCOR.err, indexARMA.err, indexALL.err)

  MAE <- abs(ERRORS.index)
  MSE <- ERRORS.index^2

  # Compute MAE metric:

  TABtot_index[, l] <- colMeans(MAE)


  # Compute RMSE metric:

  TABtot2_index[, l] <- sqrt(colMeans(MSE))

  # Take function from the empirical density of Y:

  Y0t <- Ytot.all[(train.l + 1 + h[l]):nrow(Ytot.all), 1]
  d <- density(Y0t)
  density_fun <- approxfun(d$x, d$y, method = "constant")
  w <- density_fun(Y0t)
  w <- (w - min(w)) / (max(w) - min(w))

  # Inference with absolute errors:

  TABpval[[l]] <- matrix(NA, nrow = 6, ncol = 4)

  for (i in 1:nrow(TABpval[[l]])) {
    derr <- MAE[, 1] - MAE[, 4 + i]
    derr <- w * derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval[[l]][i, 1] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)

    derr <- MAE[, 2] - MAE[, 4 + i]
    derr <- w * derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval[[l]][i, 2] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)

    derr <- MAE[, 3] - MAE[, 4 + i]
    derr <- w * derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval[[l]][i, 3] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)

    derr <- MAE[, 4] - MAE[, 4 + i]
    derr <- w * derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval[[l]][i, 4] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)
  }

  colnames(TABpval[[l]]) <- c("Base", "RW", "BU", "NO CLUST")
  rownames(TABpval[[l]]) <- c("IND", "MRKT", "EUCL", "COR", "ARMA", "ALL")

  # Inference with squared errors:

  TABpval2[[l]] <- matrix(NA, nrow = 6, ncol = 4)

  for (i in 1:nrow(TABpval2[[l]])) {
    derr <- MSE[, 1] - MSE[, 4 + i]
    derr <- w * derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval2[[l]][i, 1] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)

    derr <- MSE[, 2] - MSE[, 4 + i]
    derr <- w * derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval2[[l]][i, 2] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)

    derr <- MSE[, 3] - MSE[, 4 + i]
    derr <- w * derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval2[[l]][i, 3] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)

    derr <- MSE[, 4] - MSE[, 4 + i]
    derr <- w * derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval2[[l]][i, 4] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)
  }

  colnames(TABpval2[[l]]) <- c("Base", "RW", "BU", "NO CLUST")
  rownames(TABpval2[[l]]) <- c("IND", "MRKT", "EUCL", "COR", "ARMA", "ALL")

  ################## COSTITUENTS FORECASTING (BOTTOM SERIES)

  # 1) BASE FORECASTS:

  # Are the same of bottom-up:

  fbase.err <- matrix(NA, nrow = length(IND.f), ncol = 30)
  fbase.perr <- matrix(NA, nrow = length(IND.f), ncol = 30)
  for (i in 1:nrow(fbase.err)) {
    fbase.err[i, ] <- Ytot.I[(train.l + i + h[l] - 1), 2:ncol(Ytot.I)] - BU.rec[[i]][-1, h[l]] # BU.rec -1 to remove the top series
  }
  for (i in 1:nrow(fbase.err)) {
    fbase.perr[i, ] <- 100 * (Ytot.I[(train.l + i + h[l] - 1), 2:ncol(Ytot.I)] - BU.rec[[i]][-1, h[l]]) / Ytot.I[(train.l + i + h[l] - 1), 2:ncol(Ytot.I)]
  }
  mspe.base <- apply(fbase.err, 2, function(x) {
    sqrt(mean(x^2))
  }) # log(x^2+1)
  mape.base <- apply(fbase.err, 2, function(x) {
    mean(abs(x))
  }) # log(abs(x+1))
  mean(mape.base)
  mean(mspe.base)

  # 2) RW:

  RW.err <- matrix(NA, nrow = length(IND.f), ncol = 30)
  RW.perr <- matrix(NA, nrow = length(IND.f), ncol = 30)
  for (i in 1:nrow(RW.err)) {
    RW.err[i, ] <- Ytot.I[(train.l + i + h[l] - 1), -1] - RW.f[[i]][[1]][h[l], -1]
  }
  for (i in 1:nrow(RW.err)) {
    RW.perr[i, ] <- 100 * (Ytot.I[(train.l + i + h[l] - 1), -1] - RW.f[[i]][[1]][h[l], -1]) / Ytot.I[(train.l + i + h[l] - 1), -1]
  }
  mspe.RW <- apply(RW.err, 2, function(x) {
    sqrt(mean(x^2))
  })
  mape.RW <- apply(RW.err, 2, function(x) {
    mean(abs(x))
  })
  mean(mape.RW)
  mean(mspe.RW)

  # 3) No clustering

  recNOCLUST.err <- matrix(NA, nrow = length(ALL.f), ncol = 30)
  recNOCLUST.perr <- matrix(NA, nrow = length(ALL.f), ncol = 30)
  for (i in 1:nrow(fbase.err)) {
    recNOCLUST.err[i, ] <- Ytot.I[(train.l + i + h[l] - 1), 2:ncol(Ytot.I)] - NOCLUST.rec[[i]][1, ]
  }
  for (i in 1:nrow(fbase.err)) {
    recNOCLUST.perr[i, ] <- 100 * (Ytot.I[(train.l + i + h[l] - 1), 2:ncol(Ytot.I)] - NOCLUST.rec[[i]][1, ]) / Ytot.I[(train.l + i + h[l] - 1), 2:ncol(Ytot.I)]
  }
  mspeNOCLUST.base <- apply(recNOCLUST.err, 2, function(x) {
    sqrt(mean(x^2))
  })
  mapeNOCLUST.base <- apply(recNOCLUST.err, 2, function(x) {
    mean(abs(x))
  })
  mean(mapeNOCLUST.base)
  mean(mspeNOCLUST.base)

  # 4) IND-based reconcilied forecasts:

  recIND.err <- matrix(NA, nrow = length(IND.f), ncol = 30)
  recIND.perr <- matrix(NA, nrow = length(IND.f), ncol = 30)
  for (i in 1:nrow(fbase.err)) {
    recIND.err[i, ] <- Ytot.IND[(train.l + i + h[l] - 1), (nclust.IND + 2):ncol(Ytot.IND)] - IND.rec[[i]][1, ]
  }
  for (i in 1:nrow(fbase.err)) {
    recIND.perr[i, ] <- 100 * (Ytot.IND[(train.l + i + h[l] - 1), (nclust.IND + 2):ncol(Ytot.IND)] - IND.rec[[i]][1, ]) / Ytot.IND[(train.l + i + h[l] - 1), (nclust.IND + 2):ncol(Ytot.IND)]
  }
  mspeIND.base <- apply(recIND.err, 2, function(x) {
    sqrt(mean(x^2))
  })
  mapeIND.base <- apply(recIND.err, 2, function(x) {
    mean(abs(x))
  })
  mean(mapeIND.base)
  mean(mspeIND.base)

  # 5) MRKT-based reconcilied forecasts:

  recMRKT.err <- matrix(NA, nrow = length(MRKT.f), ncol = 30)
  recMRKT.perr <- matrix(NA, nrow = length(MRKT.f), ncol = 30)
  for (i in 1:nrow(fbase.err)) {
    recMRKT.err[i, ] <- Ytot.MRKT[(train.l + i + h[l] - 1), (nclust.MRKT + 2):ncol(Ytot.MRKT)] - MRKT.rec[[i]][1, ]
  }
  for (i in 1:nrow(fbase.err)) {
    recMRKT.perr[i, ] <- 100 * (Ytot.MRKT[(train.l + i + h[l] - 1), (nclust.MRKT + 2):ncol(Ytot.MRKT)] - MRKT.rec[[i]][1, ]) / Ytot.MRKT[(train.l + i + h[l] - 1), (nclust.MRKT + 2):ncol(Ytot.MRKT)]
  }
  mspeMRKT.base <- apply(recMRKT.err, 2, function(x) {
    sqrt(mean(x^2))
  })
  mapeMRKT.base <- apply(recMRKT.err, 2, function(x) {
    mean(abs(x))
  })
  mean(mapeMRKT.base)
  mean(mspeMRKT.base)

  # 6) EUCL-based reconcilied forecasts:

  recEUCL.err <- matrix(NA, nrow = length(EUCL.f), ncol = 30)
  recEUCL.perr <- matrix(NA, nrow = length(EUCL.f), ncol = 30)
  for (i in 1:nrow(fbase.err)) {
    recEUCL.err[i, ] <- Ytot.EUCL[(train.l + i + h[l] - 1), (nclust.EUCL + 2):ncol(Ytot.EUCL)] - EUCL.rec[[i]][1, ]
  }
  for (i in 1:nrow(fbase.err)) {
    recEUCL.perr[i, ] <- 100 * (Ytot.EUCL[(train.l + i + h[l] - 1), (nclust.EUCL + 2):ncol(Ytot.EUCL)] - EUCL.rec[[i]][1, ]) / Ytot.EUCL[(train.l + i + h[l] - 1), (nclust.EUCL + 2):ncol(Ytot.EUCL)]
  }
  mspeEUCL.base <- apply(recEUCL.err, 2, function(x) {
    sqrt(mean(x^2))
  })
  mapeEUCL.base <- apply(recEUCL.err, 2, function(x) {
    mean(abs(x))
  })
  mean(mapeEUCL.base)
  mean(mspeEUCL.base)

  # 7) COR-based reconcilied forecasts:

  recCOR.err <- matrix(NA, nrow = length(COR.f), ncol = 30)
  recCOR.perr <- matrix(NA, nrow = length(COR.f), ncol = 30)
  for (i in 1:nrow(fbase.err)) {
    recCOR.err[i, ] <- Ytot.COR[(train.l + i + h[l] - 1), (nclust.COR + 2):ncol(Ytot.COR)] - COR.rec[[i]][1, ]
  }
  for (i in 1:nrow(fbase.err)) {
    recCOR.perr[i, ] <- 100 * (Ytot.COR[(train.l + i + h[l] - 1), (nclust.COR + 2):ncol(Ytot.COR)] - COR.rec[[i]][1, ]) / Ytot.COR[(train.l + i + h[l] - 1), (nclust.COR + 2):ncol(Ytot.COR)]
  }
  mspeCOR.base <- apply(recCOR.err, 2, function(x) {
    sqrt(mean(x^2))
  })
  mapeCOR.base <- apply(recCOR.err, 2, function(x) {
    mean(abs(x))
  })
  mean(mapeCOR.base)
  mean(mspeCOR.base)

  # 8) ARMA-based reconcilied forecasts:

  orderARMA <- c(
    "AAPL", "AMGN", "DOW", "JNJ", "KO", "MRK", "MSFT", "PG", "CRM", "CSCO", "HON", "IBM", "MMM", "GS", "HD",
    "MCD", "WMT", "TRV", "VZ", "AXP", "BA", "CAT", "CVX", "DIS", "INTC", "JPM", "NKE", "UNH", "V", "WBA"
  )

  recARMA.err <- matrix(NA, nrow = length(ARMA.f), ncol = 30)
  recARMA.perr <- matrix(NA, nrow = length(ARMA.f), ncol = 30)
  for (i in 1:nrow(fbase.err)) {
    recARMA.err[i, ] <- Ytot.ARMA[(train.l + i + h[l] - 1), orderARMA] - ARMA.rec[[i]][1, ]
  }
  for (i in 1:nrow(fbase.err)) {
    recARMA.perr[i, ] <- 100 * (Ytot.ARMA[(train.l + i + h[l] - 1), orderARMA] - ARMA.rec[[i]][1, ]) / Ytot.ARMA[(train.l + i + h[l] - 1), orderARMA]
  }
  mspeARMA.base <- apply(recARMA.err, 2, function(x) {
    sqrt(mean(x^2))
  })
  mapeARMA.base <- apply(recARMA.err, 2, function(x) {
    mean(abs(x))
  })
  mean(mapeARMA.base)
  mean(mspeARMA.base)

  # 9) ALL-based reconcilied forecasts:

  recALL.err <- matrix(NA, nrow = length(ALL.f), ncol = 30)
  recALL.perr <- matrix(NA, nrow = length(ALL.f), ncol = 30)
  for (i in 1:nrow(fbase.err)) {
    recALL.err[i, ] <- ifelse(h[l] == 1, Ytot.all[(train.l + i + h[l] - 1), (nclust.all + 2):ncol(Ytot.all)] - ALL.rec[[i]], Ytot.all[(train.l + i + h[l] - 1), (nclust.all + 2):ncol(Ytot.all)] - ALL.rec[[i]][, h[l]]) # ALL.rec[[i]] if h=1, ALL.rec[[i]][,h[l]] for h>1.
  }
  for (i in 1:nrow(fbase.err)) {
    recALL.perr[i, ] <- ifelse(h[l] == 1, 100 * (Ytot.all[(train.l + i + h[l] - 1), (nclust.all + 2):ncol(Ytot.all)] - ALL.rec[[i]]) / Ytot.all[(train.l + i + h[l] - 1), (nclust.all + 2):ncol(Ytot.all)], 100 * (Ytot.all[(train.l + i + h[l] - 1), (nclust.all + 2):ncol(Ytot.all)] - ALL.rec[[i]][, h[l]]) / Ytot.all[(train.l + i + h[l] - 1), (nclust.all + 2):ncol(Ytot.all)])
  }
  mspeALL.base <- apply(recALL.err, 2, function(x) {
    sqrt(mean(x^2))
  })
  mapeALL.base <- apply(recALL.err, 2, function(x) {
    mean(abs(x))
  })
  mean(mapeALL.base)
  mean(mspeALL.base)

  # Accuracy measures comparison - average

  # MAE:

  MAE <- cbind(
    rowMeans(abs(fbase.perr)), rowMeans(abs(RW.perr)), rowMeans(abs(recNOCLUST.perr)),
    rowMeans(abs(recIND.perr)), rowMeans(abs(recMRKT.perr)), rowMeans(abs(recEUCL.perr)),
    rowMeans(abs(recCOR.perr)), rowMeans(abs(recARMA.perr)), rowMeans(abs(recALL.perr))
  )

  MAE.BASE <- rbind(
    mean(mape.base), mean(mape.RW), mean(mapeNOCLUST.base),
    mean(mapeIND.base), mean(mapeMRKT.base),
    mean(mapeEUCL.base), mean(mapeCOR.base),
    mean(mapeARMA.base), mean(mapeALL.base)
  )

  names(MAE.BASE) <- c("Base", "RW", "NO CLUST", "IND", "MRKT", "EUCL", "COR", "ARMA", "ALL")

  TABtot_bottom[, l] <- MAE.BASE

  # MSE:

  MSE <- cbind(
    rowMeans(fbase.perr^2), rowMeans(RW.perr^2), rowMeans(recNOCLUST.perr^2),
    rowMeans(recIND.perr^2), rowMeans(recMRKT.perr^2), rowMeans(recEUCL.perr^2),
    rowMeans(recCOR.perr^2), rowMeans(recARMA.perr^2), rowMeans(recALL.perr^2)
  )

  MSE.BASE <- rbind(
    mean(mspe.base), mean(mspe.RW), mean(mspeNOCLUST.base),
    mean(mspeIND.base), mean(mspeMRKT.base),
    mean(mspeEUCL.base), mean(mspeCOR.base),
    mean(mspeARMA.base), mean(mspeALL.base)
  )

  names(MSE.BASE) <- c("Base", "RW", "NO CLUST", "IND", "MRKT", "EUCL", "COR", "ARMA", "ALL")

  TABtot2_bottom[, l] <- MSE.BASE

  # Take function from the empirical density of Y:

  Y0t <- Ytot.all[(train.l + 1 + h[l]):nrow(Ytot.all), 1]
  d <- density(Y0t)
  density_fun <- approxfun(d$x, d$y, method = "constant")
  w <- density_fun(Y0t)
  w <- (w - min(w)) / (max(w))

  # Inference with MAE loss:

  TABpval.bottom[[l]] <- matrix(NA, nrow = 6, ncol = 3)

  for (i in 1:nrow(TABpval.bottom[[l]])) {
    derr <- MAE[, 1] - MAE[, 3 + i]
    # derr <- w*derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval.bottom[[l]][i, 1] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)

    derr <- MAE[, 2] - MAE[, 3 + i]
    # derr <- w*derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval.bottom[[l]][i, 2] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)

    derr <- MAE[, 3] - MAE[, 3 + i]
    # derr <- w*derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval.bottom[[l]][i, 3] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)
  }

  colnames(TABpval.bottom[[l]]) <- c("Base", "RW", "NO CLUST")
  rownames(TABpval.bottom[[l]]) <- c("IND", "MRKT", "EUCL", "COR", "ARMA", "ALL")

  # Inference with RMSE loss:

  TABpval2.bottom[[l]] <- matrix(NA, nrow = 6, ncol = 3)

  for (i in 1:nrow(TABpval2.bottom[[l]])) {
    derr <- MSE[, 1] - MSE[, 3 + i]
    # derr <- w*derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval2.bottom[[l]][i, 1] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)

    derr <- MSE[, 2] - MSE[, 3 + i]
    # derr <- w*derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval2.bottom[[l]][i, 2] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)

    derr <- MSE[, 3] - MSE[, 3 + i]
    # derr <- w*derr
    lm0 <- lm(derr ~ 1)
    derr.test <- coeftest(lm0, vcov = NeweyWest(lm0, lag = 100))
    TABpval2.bottom[[l]][i, 3] <- pt(abs(derr.test[1] / derr.test[2]), df = length(derr) - 1, lower.tail = FALSE)
  }

  colnames(TABpval2.bottom[[l]]) <- c("Base", "RW", "NO CLUST")
  rownames(TABpval2.bottom[[l]]) <- c("IND", "MRKT", "EUCL", "COR", "ARMA", "ALL")
}

TABtot_index
TABtot2_index

xtable::xtable(TABtot_index, digits = 2)
xtable::xtable(TABtot2_index, digits = 2)

TABtot_bottom
TABtot_bottom

xtable::xtable(TABtot_bottom, digits = 2)
xtable::xtable(TABtot2_bottom, digits = 2)


TABpval
TABpval2

TABpval.bottom
TABpval2.bottom

# Figures 14 and 15:

MAEbar <- matrix(NA, 6, 4)
RMSEbar <- matrix(NA, 6, 4)

for (l in 1:4) {
  for (i in 1:6) {
    MAEbar[i, l] <- log(TABtot_index[1, l] / TABtot_index[4 + i, l])
    RMSEbar[i, l] <- log(TABtot2_index[1, l] / TABtot2_index[4 + i, l])
  }
}

library(RColorBrewer)
coul <- brewer.pal(8, "Set3")

pdf(here("Figures/Best MAE2.pdf"), width = 8, height = 5)
par(mfrow = c(2, 2), mar = c(3, 3, 3, 3))
barplot(MAEbar[, 1], main = "Forecasting horizon h=1", col = coul)
barplot(MAEbar[, 2], main = "Forecasting horizon h=3", col = coul)
barplot(MAEbar[, 3], main = "Forecasting horizon h=6", col = coul)
barplot(MAEbar[, 4], main = "Forecasting horizon h=12", col = coul)
crop::dev.off.crop()

pdf(here("Figures/Best RMSE2.pdf"), width = 8, height = 5)
par(mfrow = c(2, 2), mar = c(3, 3, 3, 3))
barplot(RMSEbar[, 1], main = "Forecasting horizon h=1", col = coul)
barplot(RMSEbar[, 2], main = "Forecasting horizon h=3", col = coul)
barplot(RMSEbar[, 3], main = "Forecasting horizon h=6", col = coul)
barplot(RMSEbar[, 4], main = "Forecasting horizon h=12", col = coul)
crop::dev.off.crop()


pdf(here("Figures/Best_method.pdf"), width = 9, height = 4)
rbind(MAEbar, RMSEbar) |> 
  as.data.frame() |> 
  mutate(
    Method = rep(c("IND","MRKT","EUCL","COR","ARMA","ALL"), 2),
    Measure = rep(c("MAE","RMSE"),rep(6,2))
  ) |> 
  rename(h1 = V1, h3 = V2, h6 = V3, h12 = V4) |>
  pivot_longer(h1:h12, values_to = "Value", names_to = "Horizon") |> 
  mutate(Horizon = readr::parse_number(Horizon)) |> 
  ggplot(aes(x = Horizon, y = Value, group = Method, col = Method)) +
  geom_line() +
  facet_grid(. ~ Measure) +
  theme_bw()
crop::dev.off.crop()
