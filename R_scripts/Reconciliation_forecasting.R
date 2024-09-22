#
# Code for reproducing forecasting results in Section 3 and Section 4
#
###########################################################

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

source("Reconciliation_clustering.R")

### B: Forecasting reconciliation

train.l <- 400
out.l <- nrow(DataP0) - train.l
h <- c(1, 3, 6, 12)

Ylists <- list(Ytot.I, Ytot.IND, Ytot.MRKT, Ytot.EUCL, Ytot.COR, Ytot.ARMA)
Hlists <- list("NOCLUST" = Hy.I, "IND" = Hy.IND, "MRKT" = Hy.MRKT, "EUCL" = Hy.EUCL, "COR" = Hy.COR, "ARMA" = Hy.ARMA)

# Forecasting exercises begins:

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Set the RDS folder where forecasts will be saved:

folder <- paste0(here("RDSForecast"),"/")

# Cluster-based reconciliation: alternative approaches

for (l in 1:length(h)) {
  for (m in 1:6) {
    Ytot <- Ylists[[m]]
    Hy <- Hlists[[m]]

    for.h <- foreach(i = 1:(out.l - h[l]), .packages = c("forecast")) %dopar% {
      yf.h <- matrix(NA, h[l], ncol(Ytot))
      resi <- matrix(NA, train.l, ncol(Ytot)) # train.l+i-1 and Ytot[1:(399+i),j] for expanding window

      # ARIMA-based forecasts:

      for (j in 1:ncol(Ytot)) {
        yf.h[, j] <- forecast(auto.arima(Ytot[i:(399 + i), j]), h = h[l])$mean
        resi[, j] <- residuals(auto.arima(Ytot[i:(399 + i), j]))
      }

      return(list(yf.h, resi))
    }

    export(for.h, paste0(folder,"for.h", h[l], names(Hlists)[m], ".RDS"))

    for.rec <- foreach(i = 1:length(for.h), .packages = c("forecast", "hts")) %dopar% {
      rfor <- MinT(fcasts = for.h[[i]][[1]], nodes = get_nodes(Hy), residual = for.h[[i]][[2]], covariance = "shr", keep = "bottom")

      return(rfor)
    }

    export(for.rec, paste0(folder,"MinTrec.h", h[l], names(Hlists)[m], ".RDS"))

    print(paste0(l, "completed"))
  }
}

# ALL reconciliation:

MinTrace <- function(forvec, resimat, Smat, Cmat, keep.bottom = TRUE) {
  Wmatrix <- hts:::shrink.estim(resimat, hts:::lowerD(resimat))[[1]]

  Reconc.f <- Smat %*% solve(t(Smat) %*% solve(Wmatrix) %*% Smat) %*% t(Smat) %*% solve(Wmatrix) %*% t(forvec)

  if (keep.bottom == TRUE) {
    Reconc.b <- Reconc.f[(nrow(Cmat) + 1):nrow(Smat), ]
    names(Reconc.b) <- rownames(Smat)[(nrow(Cmat) + 1):nrow(Smat)]
    return(Reconc.b)
  } else {
    return(Reconc.f)
  }
}

for (l in 1:length(h)) {
  Ytot0 <- Ytot.all
  Cmat0 <- Cmat
  Smat0 <- Smat

  for.h <- foreach(i = 1:(out.l - h[l]), .packages = c("forecast")) %dopar% {
    yf.h <- matrix(NA, h[l], ncol(Ytot0))
    resi <- matrix(NA, train.l, ncol(Ytot0)) # train.l+i-1 and Ytot0[1:(399+i),j] for expanding window

    # ARIMA-based forecasts:

    for (j in 1:ncol(Ytot0)) {
      yf.h[, j] <- forecast(auto.arima(Ytot0[i:(399 + i), j]), h = h[l])$mean
      resi[, j] <- residuals(auto.arima(Ytot0[i:(399 + i), j]))
    }

    return(list(yf.h, resi))
  }

  export(for.h, paste0(folder,"for.h", h[l], "ALL", ".RDS"))

  for.rec <- foreach(i = 1:length(for.h), .packages = c("forecast", "hts")) %dopar% {
    rfor <- MinTrace(for.h[[i]][[1]], for.h[[i]][[2]], Smat0, Cmat0, keep.bottom = TRUE)

    return(rfor)
  }

  export(for.rec, paste0(folder,"MinTrec.h", h[l], "ALL", ".RDS"))

  print(paste0(l, "completed"))
}

# Bottom-up:

h <- c(1, 3, 6, 12)

for (l in 1:4) {
  Base.f <- import(paste0("for.h", h[l], "NOCLUST.RDS"))

  Sm <- smatrix(Hy.I)
  G <- cbind(rep(0, ncol(Sm)), diag(1, nrow = ncol(Sm), ncol = ncol(Sm)))

  for.rec <- foreach(i = 1:length(Base.f), .packages = c("forecast", "hts")) %dopar% {
    rfor <- Sm %*% G %*% t(Base.f[[i]][[1]])

    return(rfor)
  }

  export(for.rec, paste0(folder,"rec.h", h[l], "BU", ".RDS"))
}

# Random walk forecasts:

for (l in 1:length(h)) {
  Ytot <- Ytot.I

  for.h <- foreach(i = 1:(out.l - h[l]), .packages = c("forecast")) %dopar% {
    yf.h <- matrix(NA, h[l], ncol(Ytot))
    resi <- matrix(NA, train.l - 1, ncol(Ytot)) # Use this for RW

    # RW-based forecasts:

    for (j in 1:ncol(Ytot)) {
      yf.h[, j] <- rwf(Ytot[i:(399 + i), j], drift = T, h = h[l])$mean
      resi[, j] <- na.omit(resid(rwf(Ytot[i:(399 + i), j], drift = T)))
    }

    return(list(yf.h, resi))
  }

  export(for.h, paste0(folder,"forRW.h", h[l], ".RDS"))

  print(paste0(l, "completed"))
}

Ytotlist <- list(
  "NOCLUST" = Ytot.I, "IND" = Ytot.IND, "MRKT" = Ytot.MRKT, "EUCL" = Ytot.EUCL, "COR" = Ytot.COR,
  "ARMA" = Ytot.ARMA, "ALL" = Ytot.all
)

Hylist <- list(
  "NOCLUST" = Hy.I, "IND" = Hy.IND, "MRKT" = Hy.MRKT, "EUCL" = Hy.EUCL, "COR" = Hy.COR,
  "ARMA" = Hy.ARMA
)

export(Ytotlist, paste0(folder,"Ytotlist.RDS"))
export(Hylist, paste0(folder,"Hylist.RDS"))
