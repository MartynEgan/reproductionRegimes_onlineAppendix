######################################
# Passive PCA, first factorial plane #
######################################

# Load packages
library("dplyr")
library("psych")
library("ggplot2")
library("car")
library("GDAtools")

# Load data
ISSP5 <- readRDS("data/ISSP5_GDA")
ISSP4 <- readRDS("data/ISSP4_GDA")
mod <- readRDS("data/pca_mod")

# Change sign of PC1 and PC2 scores to reflect direction of Likert scale
ISSP5$PC1 <- ISSP5$PC1 *-1
ISSP5$PC2 <- ISSP5$PC2 *-1
ISSP4$PC1 <- ISSP4$PC1 *-1
ISSP4$PC2 <- ISSP4$PC2 *-1

# Drop China from ISSP4
ISSP4 <- ISSP4[ISSP4$c_alphan != "CN",]

# Function for extracting weighted mean of components/cos2 by category
m_point <- function(df, PC1 = "PC2", PC2 = "PC1", 
                    index, wt = "WEIGHT_c", 
                    na.rm = FALSE) {
  if(class(df[[index]]) == "factor") {
    n <- length(levels(droplevels(df[[index]])))
    group <- levels(droplevels(df[[index]]))
  } else {
    n <- length(unique(df[[index]][!is.na(df[[index]])]))
    group <- sort(unique(df[[index]][!is.na(df[[index]])]))
  }  
  x <- rep(0,n)
  y <- rep(0,n)
  for(i in seq_len(n)) { 
    x[i] <- weighted.mean(df[[PC1]][df[[index]] == group[i]], 
                          df[[wt]][df[[index]] == group[i]], 
                          na.rm = na.rm)
    y[i] <- weighted.mean(df[[PC2]][df[[index]] == group[i]], 
                          df[[wt]][df[[index]] == group[i]], 
                          na.rm = na.rm)
  }
  m <- data.frame(group = group,
                  x = x,
                  y = y)
  return(m)
}

# Find mean points for c_alphan
c_5_wt_mean <- m_point(ISSP5, index = "c_alphan")
c_4_wt_mean <- m_point(ISSP4, index = "c_alphan")

# Drop countries not found in both waves
c_4_wt_mean <- c_4_wt_mean[c_4_wt_mean$group %in% unique(c_5_wt_mean$group),]
c_5_wt_mean <- c_5_wt_mean[c_5_wt_mean$group %in% unique(c_4_wt_mean$group),]
c_alph <- c_4_wt_mean$group

## Calculate effect-vectors for LMEs and Nordic/CMEs and all
# Nordic
nor_4 <- c_4_wt_mean[c_4_wt_mean$group %in% c("DK", "FI", "NO", "SE"),]
nor_5 <- c_5_wt_mean[c_5_wt_mean$group %in% c("DK", "FI", "NO", "SE"),]
# CME
cme_4 <- c_4_wt_mean[c_4_wt_mean$group %in% c("DE", "JP", "SE"),]
cme_5 <- c_5_wt_mean[c_5_wt_mean$group %in% c("DE", "JP", "SE"),]
# LME
lme_4 <- c_4_wt_mean[c_4_wt_mean$group %in% c("US", "NZ", "AU", "GB-GBN"),]
lme_5 <- c_5_wt_mean[c_5_wt_mean$group %in% c("US", "NZ", "AU", "GB-GBN"),]

## Plot LME, CME and Nordic effect-vectors
png(filename = "plots/supp_pca_wt_vect_CME_LME_Nord.png",
    width = 1600, height = 600)
par(pty = "s",
    mfrow = c(1,3),
    mgp = c(1.5,0.6,0),
    mai = c(1, 1, 1, 0.5),
    cex = 2,
    cex.main = 1,
    cex.lab = 1,
    cex.axis = 1)
# plot LME plus LME mean plus mean
plot(main = "LMEs",
     NULL,
     xlab = "PC2",
     ylab = "PC1",
     xlim = c(-0.35,0.35),
     ylim = c(-0.3,0.3),
     asp = 1)
# LME mean
arrows(0,0,mean(lme_5$x-lme_4$x),mean(lme_5$y-lme_4$y),
       length = 0.2,
       code = 2,
       lty = 1,
       lwd = 2,
       col = "black")
text(mean(lme_5$x-lme_4$x),mean(lme_5$y-lme_4$y),
     pos = 2,
     label = "LME",
     cex = 0.7)
# LMEs
arrows(0,0,lme_5$x-lme_4$x,lme_5$y-lme_4$y,
       length = 0.1,
       code = 2,
       lty = 1,
       col = "grey30")
text(lme_5$x-lme_4$x,lme_5$y-lme_4$y,
     pos = c(1,2,3,2),
     label = lme_5$group,
     cex = 0.7)
# overall mean
arrows(0,0,mean(c_5_wt_mean$x-c_4_wt_mean$x),mean(c_5_wt_mean$y-c_4_wt_mean$y),
       length = 0.2,
       code = 2,
       lty = 1,
       lwd = 2,
       col = "grey80")
text(mean(c_5_wt_mean$x-c_4_wt_mean$x),mean(c_5_wt_mean$y-c_4_wt_mean$y),
     pos = 3,
     label = "Mean",
     cex = 0.7)
# plot CME mean plus CME plus overall mean
plot(main = "CMEs",
     NULL,
     xlab = "PC2",
     ylab = "PC1",
     xlim = c(-0.35,0.35),
     ylim = c(-0.3,0.3),
     asp = 1)
# CME mean
arrows(0,0,mean(cme_5$x-cme_4$x),mean(cme_5$y-cme_4$y),
       length = 0.2,
       code = 2,
       lty = 1,
       lwd = 2,
       col = "black")
text(mean(cme_5$x-cme_4$x)-0.03,mean(cme_5$y-cme_4$y)-0.03,
     label = "CME",
     cex = 0.7)
# CMEs
arrows(0,0,cme_5$x-cme_4$x,cme_5$y-cme_4$y,
       length = 0.1,
       code = 2,
       lty = 1,
       col = "grey30")
text(cme_5$x-cme_4$x,cme_5$y-cme_4$y,
     pos = c(2,2,1),
     label = cme_5$group,
     cex = 0.7)
# overall mean
arrows(0,0,mean(c_5_wt_mean$x-c_4_wt_mean$x),mean(c_5_wt_mean$y-c_4_wt_mean$y),
       length = 0.2,
       code = 2,
       lty = 1,
       lwd = 2,
       col = "grey80")
text(mean(c_5_wt_mean$x-c_4_wt_mean$x),mean(c_5_wt_mean$y-c_4_wt_mean$y),
     pos = 3,
     label = "Mean",
     cex = 0.7)
# plot Nordic mean plus Nordics plus overall mean
plot(main = "Nordic Countries",
     NULL,
     xlab = "PC2",
     ylab = "PC1",
     xlim = c(-0.35,0.35),
     ylim = c(-0.3,0.3),
     asp = 1)
# Nordic mean
arrows(0,0,mean(nor_5$x-nor_4$x),mean(nor_5$y-nor_4$y),
       length = 0.2,
       code = 2,
       lty = 1,
       lwd = 2,
       col = "black")
text(mean(nor_5$x-nor_4$x),mean(nor_5$y-nor_4$y),
     pos = 4,
     label = "Nordic",
     cex = 0.7)
# Nordics
arrows(0,0,nor_5$x-nor_4$x,nor_5$y-nor_4$y,
       length = 0.1,
       code = 2,
       lty = 1,
       col = "grey30")
text(nor_5$x-nor_4$x, 
     nor_5$y-nor_4$y,
     pos = c(2,3,4,1),
     label = nor_5$group,
     cex = 0.7)
# overall mean
arrows(0,0,mean(c_5_wt_mean$x-c_4_wt_mean$x),mean(c_5_wt_mean$y-c_4_wt_mean$y),
       length = 0.2,
       code = 2,
       lty = 1,
       lwd = 2,
       col = "grey80")
text(mean(c_5_wt_mean$x-c_4_wt_mean$x),mean(c_5_wt_mean$y-c_4_wt_mean$y),
     pos = 3,
     label = "Mean",
     cex = 0.7)
dev.off()

## Inference
# Calculate coordinates of effect points
e_points <- matrix(c(c_5_wt_mean$x-c_4_wt_mean$x, 
                     c_5_wt_mean$y-c_4_wt_mean$y),
                   nrow = nrow(c_5_wt_mean),
                   ncol = 2)
row.names(e_points) <- c_5_wt_mean$group
# Repeat for LME
lme_points <- matrix(c(lme_5$x-lme_4$x,
                       lme_5$y-lme_4$y),
                     nrow = nrow(lme_5),
                     ncol = 2)
row.names(lme_points) <- lme_5$group
# Repeat for CME
cme_points <- matrix(c(cme_5$x-cme_4$x,
                       cme_5$y-cme_4$y),
                     nrow = nrow(cme_5),
                     ncol = 2)
row.names(cme_points) <- cme_5$group
# Repeat for Nordic
nordic_points <- matrix(c(nor_5$x-nor_4$x,
                          nor_5$y-nor_4$y),
                        nrow = nrow(nor_5),
                        ncol = 2)
row.names(nordic_points) <- nor_5$group

## Perform geometric typicality test for mean effect-vector using  
## Geometric_Typicality.R script.
base <- e_points
evmcov <- Mcov.KK <- cov.wt(base, method = "ML")$cov
source("code/gt_test.R")
ev_gt <- gt.test(base)
# M distance: 0.748 (notable)
# P value: 2458/65536 = 0.038
# Kappa: 0.715 (95% confidence region between 0.698 and 0.73)

base <- ISSP4[,c("PC2","PC1")]
Mcov.KK <- cov.wt(base, wt = ISSP4$WEIGHT_c, method = "ML")$cov
ev_full_gt <- gt.test(base, centre = FALSE, mean = c(
  weighted.mean(ISSP5$PC2, w = ISSP5$WEIGHT_c),
  weighted.mean(ISSP5$PC1, w = ISSP5$WEIGHT_c)
))

## Plot effect-vectors with compatibility ellipse (geometric typicality test)
# label adjustment
pos <- c(1,2,3,3,1,1,3,2,3,3,3,1,2,3,3,4) 
png(filename = "plots/c_alphan_wt_vectors.png",
    width = 1400, height = 1400)
par(pty = "s",
    cex = 2.5)
plot(main = NULL,#"Effect-Vectors by Country, 2009-2019 \n (Supplementary PCA, First Factorial Plane)",
     NULL,
     xlab = paste0("PC2, ", round(mod$Vaccounted[[2,2]]*100, 1), "% of variance\n (Capital Structure/Symbolic System)"),
     ylab = paste0("PC1, ", round(mod$Vaccounted[[2,1]]*100, 1), "% of variance\n (Capital Volume)"),
     xlim = c(-0.35,0.35),
     ylim = c(-0.35,0.35),
     asp = 1)
arrows(0,0,c_5_wt_mean$x-c_4_wt_mean$x, c_5_wt_mean$y-c_4_wt_mean$y,
       length = 0.2,
       code = 2,
       lty = 1,
       col = "grey30")
text(c_5_wt_mean$x-c_4_wt_mean$x, c_5_wt_mean$y-c_4_wt_mean$y,
     pos = pos,
     label = c_5_wt_mean$group,
     cex = 0.8)
arrows(0,0,mean(c_5_wt_mean$x-c_4_wt_mean$x),mean(c_5_wt_mean$y-c_4_wt_mean$y),
       length = 0.2,
       code = 2,
       lty = 1,
       lwd = 4,
       col = "black")
text(mean(c_5_wt_mean$x-c_4_wt_mean$x),mean(c_5_wt_mean$y-c_4_wt_mean$y),
     pos = 2,
     label = substitute(paste(bold("MEAN"))),
     cex = 0.8)
ellipse(c(mean(c_5_wt_mean$x-c_4_wt_mean$x), 
          mean(c_5_wt_mean$y-c_4_wt_mean$y)),
        evmcov, #generate from geometric typicality test
        ev_gt$kappa.D, # kappa from geometric typicality test
        col = "grey50",
        lty = 3,
        center.pch = FALSE)
dev.off()

## Combinatorial typicality test on c_alphan clouds

# Parameters
notable_D  <- 0.34
alpha    <- 0.05
max_number <- 100000
seed       <- 2022

## Descriptive analysis (M-Distance)
# Each c_alphan
m_dist <- matrix(rep(0,length(c_alph)), ncol = 1)
row.names(m_dist) <- c_alph

for (i in 1:length(c_alph)) {
  base <- ISSP5[ISSP5$c_alphan %in% c_alph[i],c("PC2","PC1")]
  wt <- ISSP5[ISSP5$c_alphan %in% c_alph[i],"WEIGHT_c"]
  group <- ISSP4[ISSP4$c_alphan %in% c_alph[i],c("PC2","PC1")]
  wt4 <- ISSP4[ISSP4$c_alphan %in% c_alph[i],"WEIGHT_c"]
  # code from B. Le Roux
  n    <- dim(base)[1]
  K    <- dim(base)[2]
  X.IK <- as.matrix(base, nrow= n, ncol= K)
  n_c  <- dim(group)[1]
  #X.CK <- as.matrix(group, nrow= n_c, ncol= K)
  if(dim(group)[2] < K)
    stop("the number of columns for group is less than ", K)
  Mcov.KK <- cov.wt(X.IK, wt = wt, method= "ML")$cov
  eig <- eigen(Mcov.KK, symmetric = TRUE)
  L <- sum(eig$values > 1.5e-8)
  lambda.L <- eig$values[1:L]
  if (L < K) {
    warning("The dimension of the reference cloud is ", L, " < ", K,
            " (number of variables). \n It is advisable to study again",
            " the geometric construction of the reference cloud.")
  }
  BasisChange.KL <- eig$vectors[ ,1:L] %*% diag(1/sqrt(lambda.L), nrow=L)
  base.m <- sapply(base, weighted.mean, w = wt)
  Z.IL   <- sweep(X.IK, 2, base.m, "-") %*% BasisChange.KL
  group.m <- sapply(group, weighted.mean, w = wt4)
  ZC.L   <- t(group.m - base.m) %*% BasisChange.KL
  d2_obs <- sum(ZC.L^2)
  m_dist[i,] <- sqrt(d2_obs)
}  
m_dist_o <- as.matrix(m_dist[order(m_dist, decreasing = TRUE),])


## Plot Effect Vectors in Factorial Plane
# label adjustments
pos_4 <- c(4,4,4,1,3,4,1,4,1,4,1,3,4,1,3,3)
pos_5 <- c(1,2,2,3,1,2,3,3,3,3,3,1,3,3,3,3)

# m_distance circle
radius = 0.4
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) 

# Plot
png(filename = "plots/c_alphan_wt_suppPCA.png",
    width = 1200, height = 1200)
par(pty = "s", cex = 2)
plot(c_5_wt_mean$x, c_5_wt_mean$y,
     pch = 0,
     col = "black",
     xlab = paste0("PC2, ", round(mod$Vaccounted[[2,2]]*100, 1), "% of variance\n (Capital Structure/Symbolic System)"),
     ylab = paste0("PC1, ", round(mod$Vaccounted[[2,1]]*100, 1), "% of variance\n (Capital Volume)"),
     main = NULL, #"Country Transformations in First Factorial Plane \n (Supplementary PCA)",
     xlim = c(-0.8,0.8),
     ylim = c(-0.8,0.8),
     asp = 1)
text(c(0.88,0.88,-0.88,-0.88,0,0),c(0.03,-0.03,0.03,-0.03,0.88,-0.88),
     pos = c(2,2,4,4,3,1),
     label = c("MERIT", "CULTURAL CAPITAL+",
               "DISENCHANTMENT", "CULTURAL CAPITAL-",
               "CAPITAL VOLUME+", "CAPITAL VOLUME-"),
     cex = 0.8)
points(c_4_wt_mean$x, c_4_wt_mean$y,
       pch = 15,
       col = "grey70")
text(c_5_wt_mean$x, c_5_wt_mean$y,
     pos = pos_5,
     labels = ifelse(m_dist >= 0.4, paste0(round(m_dist,2),"*"),
                     ifelse(m_dist < 0.4 & m_dist >= 0.34, paste0(round(m_dist,2),"\U207A"),
                                                                  round(m_dist,2))))
text(c_4_wt_mean$x, c_4_wt_mean$y,
     pos = pos_4,
     labels = c_4_wt_mean$group)
arrows(c_4_wt_mean$x, c_4_wt_mean$y,
       c_5_wt_mean$x, c_5_wt_mean$y,
       length = 0.2,
       code = 2,
       lty = 1,
       col = "grey30")
legend("topright",
       legend = c("2009", "2019"),
       col = c("grey70", "black"), pch = c(15, 0),
       cex = 1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y,
      lty = 1, col = "grey80")
dev.off()