#
# Code for reproducing clustering results in Section 3 and Section 4
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

# Default colours
options(
  ggplot2.discrete.colour = c("#D55E00", "#0072B2","#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442"),
  ggplot2.discrete.fill = c("#D55E00", "#0072B2","#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442")
)

# Load data:

DataP0 <- import(here("DataP0.xlsx"))
DJIA.dat <- import(here("DJIA - composition at 2020.xlsx"))
Basey <- DataP0[, -1]

# Tidy version of data for plotting
DJI_df <- DataP0 |>
  mutate(
    Date = as.Date(Dates0),
    DJI = rowSums(Basey)
  ) |>
  select(-Dates0) |>
  tidyr::pivot_longer(-Date, names_to = "Stock", values_to = "Price")

# Plot time series:

pdf(here("Figures/DJIA_constitutes.pdf"), width = 8, height = 6)
#par(mfrow = c(6, 5), mar = c(2, 1, 2, 1))
#plot(DataP0[, 1], Basey[, 1], type = "l", main = colnames(Basey)[1])
#for (i in 2:ncol(Basey)) {
#  plot(DataP0[, 1], Basey[, i], type = "l", main = colnames(Basey)[i])
#}
DJI_df |>
  filter(Stock != "DJI") |>
  ggplot(aes(x = Date, y = Price)) +
  geom_line() +
  facet_wrap(facets = vars(Stock), scales = "free_y", ncol=5) +
  scale_x_date(date_breaks = "1 year", date_labels="%Y") +
  theme_bw()
crop::dev.off.crop()

pdf(here("Figures/DJIAts.pdf"), width = 8, height = 4)
#plot(DataP0[, 1], rowSums(Basey), type = "l",
#     main = "Dow Jones Industrial Average", ylab = "Price", xlab = "Year")
p1 <- DJI_df |>
  filter(Stock == "DJI") |>
  ggplot(aes(x = Date, y = Price)) +
  geom_line() +
  labs(title = "Dow Jones Industrial Average") +
  scale_x_date(date_breaks = "1 year", minor_breaks = "1 month", date_labels = "%Y") +
  theme_bw()
print(p1)
crop::dev.off.crop()

# A. Discovering clustering structure:

# 0. No clustering:

Hy.I <- hts(Basey)
Ytot.I <- aggts(Hy.I)

pdf(here("Figures/Hy_NOCLUST.pdf"), width = 8, height = 5)
#plot(Hy.I)
stocks <- sort(unique(DJI_df$Stock))
stocks <- c(stocks[stocks != "DJI"], "DJI")
DJI_df |>
  mutate(Index = if_else(Stock == "DJI", "DJI", "Stocks")) |>
  ggplot(aes(x = Date, y = Price, group = Stock, col = Stock)) +
  geom_line() +
  facet_grid(Index ~ ., scales = "free_y") +
  labs(title = "Dow Jones Industrial Average with Component Stocks") +
  scale_x_date(date_breaks = "1 year", minor_breaks = "1 month", date_labels = "%Y") +
  theme_bw() +
  scale_color_manual(breaks = stocks, values = c(scales::hue_pal()(length(stocks)-1),"#000000")) +
  guides(col = "none")
crop::dev.off.crop()

# 1. Subjective clusters:

# Industry-based clusters:

nclust.IND <- length(unique(DJIA.dat$Industry))
DJIA.dat2 <- DJIA.dat[order(DJIA.dat$Symbol), ]
clust.IND <- as.numeric(as.factor(DJIA.dat2$Industry))
DJIA.dat2 <- cbind(DJIA.dat2, clust.IND)
Hy.IND <- hts(Basey[, order(DJIA.dat2$clust.IND)], nodes = list(nclust.IND, table(clust.IND)))
Ytot.IND <- aggts(Hy.IND)

DJIA.dat2$Label = LETTERS[DJIA.dat2$clust.IND]
industry_df <- Ytot.IND |>
  as.data.frame() |>
  select(A:T) |>
  mutate(Date = as.Date(DataP0$Date)) |>
  tidyr::pivot_longer(-Date, names_to = "Label", values_to = "Price") |>
  left_join(
    DJIA.dat2 |> select(Label, Industry) |> distinct(),
    by = "Label"
  )
DJI_df <- DJI_df |>
  left_join(
    DJIA.dat2 |> select(Symbol, Label, Industry),
    by = join_by(Stock == Symbol)
  )

pdf(here("Figures/Hy_IND.pdf"), width = 8, height = 8)
#plot(Hy.IND)
industries <- sort(unique(DJIA.dat2$Industry))
bind_rows(
    DJI_df |>
      select(-Label) |>
      mutate(Group = if_else(Stock == "DJI", "DJI", "Stocks")) |>
      mutate(Industry = if_else(Stock == "DJI", "All", Industry)),
    industry_df |>
      rename(Stock = Label) |>
      mutate(Group = "Industry")
  ) |>
  mutate(Industry = factor(Industry, levels = c(industries, "All"))) |>
  ggplot(aes(x = Date, y = Price, group = Stock, col = Industry)) +
  geom_line() +
  facet_grid(Group ~ ., scales = "free_y") +
  labs(title = "Dow Jones Industrial Average with Industry Groups") +
  scale_x_date(date_breaks = "1 year", minor_breaks = "1 month", date_labels = "%Y") +
  scale_color_manual(breaks = c(industries, "All"), values = c(scales::hue_pal()(length(industries)),"#000000")) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  guides(col = guide_legend(ncol = 3)) 
crop::dev.off.crop()

# Exchange-based clusters:

nclust.MRKT <- length(unique(DJIA.dat$Exchange))
DJIA.dat2 <- DJIA.dat[order(DJIA.dat$Symbol), ]
clust.MRKT <- as.numeric(as.factor(DJIA.dat2$Exchange))
DJIA.dat2 <- cbind(DJIA.dat2, clust.MRKT)
Hy.MRKT <- hts(Basey[, order(DJIA.dat2$clust.MRKT)], nodes = list(nclust.MRKT, table(clust.MRKT)))
Ytot.MRKT <- aggts(Hy.MRKT)

DJIA.dat2$Label = LETTERS[DJIA.dat2$clust.MRKT]
market_df <- Ytot.MRKT |>
  as.data.frame() |>
  select(A:B) |>
  mutate(Date = as.Date(DataP0$Date)) |>
  tidyr::pivot_longer(-Date, names_to = "Label", values_to = "Price") |>
  left_join(
    DJIA.dat2 |> select(Label, Exchange) |> distinct(),
    by = "Label"
  )
DJI_df <- DJI_df |>
  select(-Label) |>
  left_join(
    DJIA.dat2 |> select(Symbol, Label, Exchange),
    by = join_by(Stock == Symbol)
  )

pdf(here("Figures/Hy_MRKT.pdf"), width = 8, height = 7)
#plot(Hy.MRKT)
exchanges <- sort(unique(DJIA.dat2$Exchange))
bind_rows(
    DJI_df |>
      select(-Label) |>
      mutate(Group = if_else(Stock == "DJI", "DJI", "Stocks")) |>
      mutate(Exchange = if_else(Stock == "DJI", "All", Exchange)),
    market_df |>
      rename(Stock = Label) |>
      mutate(Group = "Exchange")
  ) |>
  mutate(Exchange = factor(Exchange, levels = c(exchanges, "All"))) |>
  ggplot(aes(x = Date, y = Price, group = Stock, col = Exchange)) +
  geom_line() +
  facet_grid(Group ~ ., scales = "free_y") +
  labs(title = "Dow Jones Industrial Average with Exchange Groups") +
  scale_x_date(date_breaks = "1 year", minor_breaks = "1 month", date_labels = "%Y") +
  scale_color_manual(breaks = c(exchanges, "All"), values = c("#D55E00", "#0072B2","#000000")) +
  theme_bw() +
  theme(legend.position = "bottom")
crop::dev.off.crop()


# 2. Unsupervised approaches:

train.l <- 400

# Returns-based measures:

# EUCL-based clusters (on returns):

DataP0 <- import(here("DataP0.xlsx"))
Basey <- DataP0[, -1]
Baseret <- apply(Basey, 2, function(x) {
  diff(log(x))
})

ASW <- NULL
for (c in 1:9) {
  ASW[c] <- pam(diss(t(Baseret[1:train.l, ]), "EUCL"), k = c + 1, diss = TRUE)$silinfo$avg.width
}
pdf(here("Figures/ASW_EUCL.pdf"), width = 8, height = 3)
#plot(2:10, ASW, type = "b", main = "ASW for EUCL-based PAM clustering", ylab = "ASW", xlab = "Clusters")
data.frame(Clusters = 2:10,
           ASW = ASW) |>
  ggplot(aes(x=Clusters, y = ASW)) +
  geom_line() + geom_point() +
  labs(title = "ASW for EUCL-based PAM clustering") +
  theme_bw()
crop::dev.off.crop()

nclust.EUCL <- which.max(ASW) + 1
clust.EUCL <- pam(diss(t(Baseret[1:train.l, ]), "EUCL"), k = nclust.EUCL, diss = TRUE)
Basey.EUCL <- Basey[, order(clust.EUCL$clustering)]
Hy.EUCL <- hts(Basey.EUCL, nodes = list(nclust.EUCL, table(clust.EUCL$clustering)))
Ytot.EUCL <- aggts(Hy.EUCL)

eucl_df <- Ytot.EUCL |>
  as.data.frame() |>
  select(A:C) |>
  mutate(Date = as.Date(DataP0$Date)) |>
  tidyr::pivot_longer(-Date, names_to = "Label", values_to = "Price")
DJI_df <- DJI_df |>
  select(-Label) |>
  left_join(
    data.frame(Stock = names(clust.EUCL$clustering), Cluster = clust.EUCL$clustering) |>
      mutate(Label = LETTERS[Cluster]) |>
      select(Stock, Label),
    by = "Stock"
  ) |>
  rename(EUCL = Label)
pdf(here("Figures/Hy_EUCL.pdf"), width = 8, height = 7)
#plot(Hy.EUCL)
bind_rows(
  DJI_df |>
    mutate(Group = if_else(Stock == "DJI", "DJI", "Stocks")) |>
    mutate(Cluster = if_else(Stock == "DJI", "All", EUCL)),
  eucl_df |>
    rename(Stock = Label) |>
    mutate(Group = "Cluster", Cluster = Stock)
  ) |>
  mutate(
    Cluster = factor(Cluster, levels = c(LETTERS[1:3],"All")),
    Group = factor(Group, levels = c("DJI", "Cluster", "Stocks"))
  ) |>
  ggplot(aes(x = Date, y = Price, group = Stock, col = Cluster)) +
  geom_line() +
  facet_grid(Group ~ ., scales = "free_y") +
  labs(title = "Dow Jones Industrial Average with EUCL Cluster Groups") +
  scale_x_date(date_breaks = "1 year", minor_breaks = "1 month", date_labels = "%Y") +
  scale_color_manual(breaks = c(LETTERS[1:3], "All"), values = c("#D55E00", "#0072B2","#009E73","#000000")) +
  theme_bw() +
  theme(legend.position = "bottom")
crop::dev.off.crop()

# COR-based clusters (on returns):

DataP0 <- import(here("DataP0.xlsx"))
Basey <- DataP0[, -1]
Baseret <- apply(Basey, 2, function(x) {
  diff(log(x))
})

ASW <- NULL
for (c in 1:9) {
  ASW[c] <- pam(diss(t(Baseret[1:train.l, ]), "COR"), k = c + 1, diss = TRUE)$silinfo$avg.width
}

pdf(here("Figures/ASW_COR.pdf"), width = 8, height = 3)
#plot(2:10, ASW, type = "b", main = "ASW for COR-based PAM clustering", ylab = "ASW", xlab = "Clusters")
data.frame(Clusters = 2:10,
           ASW = ASW) |>
  ggplot(aes(x=Clusters, y = ASW)) +
  geom_line() + geom_point() +
  labs(title = "ASW for COR-based PAM clustering") +
  theme_bw()
crop::dev.off.crop()

nclust.COR <- which.max(ASW) + 1
clust.COR <- pam(diss(t(Baseret[1:train.l, ]), "COR"), k = nclust.COR, diss = TRUE)
Basey.COR <- Basey[, order(clust.COR$clustering)]
Hy.COR <- hts(Basey.COR, nodes = list(nclust.COR, table(clust.COR$clustering)))
Ytot.COR <- aggts(Hy.COR)

cor_df <- Ytot.COR |>
  as.data.frame() |>
  select(A:C) |>
  mutate(Date = as.Date(DataP0$Date)) |>
  tidyr::pivot_longer(-Date, names_to = "Label", values_to = "Price")
DJI_df <- DJI_df |>
  left_join(
    data.frame(Stock = names(clust.COR$clustering), Cluster = clust.COR$clustering) |>
      mutate(Label = LETTERS[Cluster]) |>
      select(Stock, Label),
    by = "Stock"
  ) |>
  rename(COR = Label)

pdf(here("Figures/Hy_COR.pdf"), width = 8, height = 7)
#plot(Hy.COR)
bind_rows(
  DJI_df |>
    mutate(Group = if_else(Stock == "DJI", "DJI", "Stocks")) |>
    mutate(Cluster = if_else(Stock == "DJI", "All", COR)),
  cor_df |>
    rename(Stock = Label) |>
    mutate(Group = "Cluster", Cluster = Stock)
) |>
  mutate(
    Cluster = factor(Cluster, levels = c(LETTERS[1:3],"All")),
    Group = factor(Group, levels = c("DJI", "Cluster", "Stocks"))
  ) |>
  ggplot(aes(x = Date, y = Price, group = Stock, col = Cluster)) +
  geom_line() +
  facet_grid(Group ~ ., scales = "free_y") +
  labs(title = "Dow Jones Industrial Average with COR Cluster Groups") +
  scale_x_date(date_breaks = "1 year", minor_breaks = "1 month", date_labels = "%Y") +
  scale_color_manual(breaks = c(LETTERS[1:3], "All"), values = c("#D55E00", "#0072B2","#009E73","#000000")) +
  theme_bw() +
  theme(legend.position = "bottom")
crop::dev.off.crop()

# Prices-based measures:

# ARIMA-based clustering (prices)

DataP0 <- import(here("DataP0.xlsx"))
Basey <- DataP0[, -1]

armodels <- list()

for (j in 1:ncol(Basey)) {
  armodels[[j]] <- auto.arima(Basey[1:train.l, j])
}

ModelsAR <- matrix(NA, nrow = ncol(Basey), ncol = 1)
for (j in 1:ncol(Basey)) {
  ModelsAR[j, 1] <- paste0("ARIMA(", armodels[[j]]$arma[1], ",", 1, ",", armodels[[j]]$arma[2], ")")
}
rownames(ModelsAR) <- colnames(Basey)
xtable::xtable(ModelsAR)

k <- 10
PIcoeff <- matrix(NA, k, ncol(Basey))
for (j in 1:ncol(PIcoeff)) {
  model0 <- armodels[[j]]
  PIcoeff[, j] <- ARMAtoAR(ar = model0$model$phi, ma = model0$model$theta, lag.max = k)
}
colnames(PIcoeff) <- colnames(Basey)

which(colSums(PIcoeff) == 0)

ASW <- NULL
for (c in 1:9) {
  ASW[c] <- pam(diss(t(PIcoeff[, -which(colSums(PIcoeff) == 0)]), "EUCL"), k = c + 1, diss = TRUE)$silinfo$avg.width
}
pdf(here("Figures/ASW_ARIMA.pdf"), width = 8, height = 3)
#plot(2:10, ASW, type = "b", main = "ASW for ARIMA-based PAM clustering", ylab = "ASW", xlab = "Clusters")
data.frame(Clusters = 2:10,
           ASW = ASW) |>
  ggplot(aes(x=Clusters, y = ASW)) +
  geom_line() + geom_point() +
  labs(title = "ASW for ARIMA-based PAM clustering") +
  theme_bw()
crop::dev.off.crop()

nclust.ARMA <- which.max(ASW) + 1
clust.ARMA <- pam(diss(t(PIcoeff[, -which(colSums(PIcoeff) == 0)]), "EUCL"), k = nclust.ARMA, diss = TRUE)

Clu5 <- rep(5, length(which(colSums(PIcoeff) == 0)))
names(Clu5) <- names(which(colSums(PIcoeff) == 0))
CluARMA <- c(clust.ARMA$clustering, Clu5)

clust.ARMA$clustering <- CluARMA[sort(names(CluARMA))]

pdf(here("Figures/ARIMA_plot.pdf"), width = 8, height = 3)
# plot(PIcoeff[, clust.ARMA$clustering == 1][, 1], ylim = c(min(PIcoeff), max(PIcoeff)), type = "l", main = "ARMA-based clustering: AR(∞) weights", ylab = "Coeff.", xlab = "K")
# for (i in 2:table(clust.ARMA$clustering)[1]) {
#   lines(PIcoeff[, clust.ARMA$clustering == 1][, i], col = "black")
# }
# for (i in 1:table(clust.ARMA$clustering)[2]) {
#   lines(PIcoeff[, clust.ARMA$clustering == 2][, i], col = "red")
# }
# for (i in 1:table(clust.ARMA$clustering)[3]) {
#   lines(PIcoeff[, clust.ARMA$clustering == 3][, i], col = "blue")
# }
# for (i in 1:table(clust.ARMA$clustering)[4]) {
#   lines(PIcoeff[, clust.ARMA$clustering == 4][, i], col = "green")
# }
# for (i in 1:table(clust.ARMA$clustering)[5]) {
#   lines(PIcoeff[, clust.ARMA$clustering == 5][, i], col = "orange")
# }
# legend("bottomright", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col = c("black", "red", "blue", "green", "orange"), lty = 1, cex = 0.6)
clusters <- data.frame(Stock = names(clust.ARMA$clustering), Cluster = LETTERS[clust.ARMA$clustering])
PIcoeff |>
  as.data.frame() |>
  mutate(lag = seq(NROW(PIcoeff))) |>
  tidyr::pivot_longer(-lag, names_to = "Stock", values_to = "Coefficient") |>
  left_join(clusters, by = "Stock") |>
  ggplot(aes(x = lag, y = Coefficient, group = Stock, col = Cluster)) +
  geom_line() +
  theme_bw() +
  labs(x = "lag (K)")
crop::dev.off.crop()

Basey.ARMA <- Basey[, order(clust.ARMA$clustering)]
Hy.ARMA <- hts(Basey.ARMA, nodes = list(5, table(clust.ARMA$clustering)))
Ytot.ARMA <- aggts(Hy.ARMA)

arma_df <- Ytot.ARMA |>
  as.data.frame() |>
  select(A:E) |>
  mutate(Date = as.Date(DataP0$Date)) |>
  tidyr::pivot_longer(-Date, names_to = "Label", values_to = "Price")
DJI_df <- DJI_df |>
  left_join(
    clusters |>
      mutate(Label = LETTERS[as.numeric(Cluster)]) |>
      select(Stock, Label),
    by = "Stock"
  ) |>
  rename(ARMA = Label)

pdf(here("Figures/Hy_ARMA.pdf"), width = 8, height = 7)
#plot(Hy.ARMA)
bind_rows(
  DJI_df |>
    mutate(Group = if_else(Stock == "DJI", "DJI", "Stocks")) |>
    mutate(Cluster = if_else(Stock == "DJI", "All", ARMA)),
  arma_df |>
    rename(Stock = Label) |>
    mutate(Group = "Cluster", Cluster = Stock)
) |>
  mutate(
    Cluster = factor(Cluster, levels = c(LETTERS[1:5],"All")),
    Group = factor(Group, levels = c("DJI", "Cluster", "Stocks"))
  ) |>
  ggplot(aes(x = Date, y = Price, group = Stock, col = Cluster)) +
  geom_line() +
  facet_grid(Group ~ ., scales = "free_y") +
  labs(title = "Dow Jones Industrial Average with ARMA Cluster Groups") +
  scale_x_date(date_breaks = "1 year", minor_breaks = "1 month", date_labels = "%Y") +
  scale_color_manual(breaks = c(LETTERS[1:5], "All"),
                     values = c("#D55E00", "#0072B2","#009E73", "#CC79A7", "#E69F00","#000000")) +
  theme_bw() +
  theme(legend.position = "bottom")
crop::dev.off.crop()

# Combination approaches: ALL

DataP0 <- import(here("DataP0.xlsx"))
DJIA.dat <- import(here("DJIA - composition at 2020.xlsx"))
DJIA.dat2 <- DJIA.dat[order(DJIA.dat$Symbol), ]
Basey <- DataP0[, -1]

nclust.IND <- length(unique(DJIA.dat2$Industry))
clust.IND <- as.numeric(as.factor(DJIA.dat2$Industry))
DJIA.dat2 <- cbind(DJIA.dat2, clust.IND)
nclust.MRKT <- length(unique(DJIA.dat2$Exchange))
clust.MRKT <- as.numeric(as.factor(DJIA.dat2$Exchange))
DJIA.dat2 <- cbind(DJIA.dat2, clust.MRKT)
Baseret <- apply(Basey, 2, function(x) {
  diff(log(x))
})
nclust.EUCL <- 3
clust.EUCL <- pam(diss(t(Baseret[1:train.l, ]), "EUCL"), k = nclust.EUCL, diss = TRUE)
DJIA.dat2 <- cbind(DJIA.dat2, clust.EUCL$clustering)
nclust.COR <- 3
clust.COR <- pam(diss(t(Baseret[1:train.l, ]), "COR"), k = nclust.COR, diss = TRUE)
DJIA.dat2 <- cbind(DJIA.dat2, clust.COR$clustering)

nclust.ARMA <- 4
clust.ARMA <- pam(diss(t(PIcoeff[, -which(colSums(PIcoeff) == 0)]), "EUCL"), k = nclust.ARMA, diss = TRUE)
Clu5 <- rep(5, length(which(colSums(PIcoeff) == 0)))
names(Clu5) <- names(which(colSums(PIcoeff) == 0))
CluARMA <- c(clust.ARMA$clustering, Clu5)
clust.ARMA$clustering <- CluARMA[sort(names(CluARMA))]
DJIA.dat2 <- cbind(DJIA.dat2, clust.ARMA$clustering)

nclust.all <- nclust.IND + nclust.MRKT + nclust.EUCL + nclust.COR + nclust.ARMA + 1

Cmat <- matrix(NA, ncol = ncol(Basey), nrow = nclust.all + 1)
colnames(Cmat) <- colnames(Basey)
Cmat[1, ] <- rep(1, ncol(Cmat))
for (c in 1:20) {
  Cmat[1 + c, ] <- ifelse(DJIA.dat2$clust.IND == c, 1, 0)
}
for (c in 1:2) {
  Cmat[21 + c, ] <- ifelse(DJIA.dat2$clust.MRKT == c, 1, 0)
}
for (c in 1:3) {
  Cmat[23 + c, ] <- ifelse(DJIA.dat2$`clust.EUCL$clustering` == c, 1, 0)
}
for (c in 1:3) {
  Cmat[26 + c, ] <- ifelse(DJIA.dat2$`clust.COR$clustering` == c, 1, 0)
}
for (c in 1:5) {
  Cmat[29 + c, ] <- ifelse(DJIA.dat2$`clust.ARMA$clustering` == c, 1, 0)
}
rownames(Cmat) <- 1:nrow(Cmat)
rownames(Cmat)[1] <- "tot"
for (c in 1:20) {
  rownames(Cmat)[1 + c] <- paste0("c", c, ".IND")
}
for (c in 1:2) {
  rownames(Cmat)[21 + c] <- paste0("c", c, ".MRKT")
}
for (c in 1:3) {
  rownames(Cmat)[23 + c] <- paste0("c", c, ".EUCL")
}
for (c in 1:3) {
  rownames(Cmat)[26 + c] <- paste0("c", c, ".COR")
}
for (c in 1:5) {
  rownames(Cmat)[29 + c] <- paste0("c", c, ".ARMA")
}

Smat <- rbind(Cmat, diag(1, ncol(Basey)))
rownames(Smat)[(nclust.all + 2):nrow(Smat)] <- colnames(Basey)
Ytot.all <- matrix(NA, nrow = nrow(Basey), ncol = nrow(Smat))
for (t in 1:nrow(Basey)) {
  Ytot.all[t, ] <- Smat %*% t(Basey[t, ])
}
Ytot.tot <- Ytot.all
