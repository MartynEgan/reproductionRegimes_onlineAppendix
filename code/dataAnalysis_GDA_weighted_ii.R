#####################################################################
# Data Analysis: GDA (weighted), second factorial plane (PCs 1 + 3) #
#####################################################################

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

source("code/ISCO08_code.r")

# Change sign of PC1 scores to reflect direction of Likert scale
ISSP5$PC1 <- ISSP5$PC1 *-1
ISSP4$PC1 <- ISSP4$PC1 *-1

# subset ISSP4 on CN
ISSP4 <- ISSP4[ISSP4$c_alphan == "CN",]

## Analysis and visualisation
# Function for plotting ellipse vertices
# 1. Using conf level
vertices <- function(x, y, level) {
  ell_info <- cov.wt(cbind(x, y))
  E <- eigen(ell_info$cov, symmetric = TRUE)  
  U <- E[[2]]  
  D <- sqrt(E[[1]])  
  dfn <- 2
  dfd <- length(x) - 1
  r <- sqrt(dfn * qf(level, dfn, dfd))
  Z <- rbind(c(r, 0), c(0, r), c(-r, 0), c(0, -r))  
  Z <- tcrossprod(Z * rep(D, each = 4), U)  
  Z <- Z + rep(ell_info$center, each = 4)
  segments(Z[1,1], Z[1,2], Z[3,1], Z[3,2], lty = 2, col = "black", lwd = 1)
  segments(Z[2,1], Z[2,2], Z[4,1], Z[4,2], lty = 2, col = "black", lwd = 1)
}

#2. Using kappa
vertices.k <- function(centre, Mcov, k, col = "black") {
  E <- eigen(Mcov, symmetric = TRUE)  
  U <- E[[2]]  
  D <- sqrt(E[[1]])  
  r <- k
  Z <- rbind(c(r, 0), c(0, r), c(-r, 0), c(0, -r))  
  Z <- tcrossprod(Z * rep(D, each = 4), U)  
  Z <- Z + rep(centre, each = 4)
  segments(Z[1,1], Z[1,2], Z[3,1], Z[3,2], lty = 2, col = col, lwd = 1)
  segments(Z[2,1], Z[2,2], Z[4,1], Z[4,2], lty = 2, col = col, lwd = 1)
}

# Function for calculating ellipse data
el.data <- function(Mcov){
  E <- eigen(Mcov, symmetric = TRUE)  
  D <- sqrt(E[[1]])  
  e <- sqrt(1-((min(D)^2)/max(D)^2))
  theta <- atan(E[[2]][[1,1]]/E[[2]][[2,1]]) * (180/pi)
  out <- list()
  out$eig <- E
  out$D <- D
  out$e <- e
  out$theta <- theta
  return(out)
}

# Function for performing geometric typicality test
source("code/gt_test.R")
by_gt <- function(df, PC1 = "PC3", PC2 = "PC1", index, 
                  weight = "WEIGHT_c", notable_D = 0.4, 
                  alpha = 0.05, max_number = 1e+06, 
                  n_dir = 200, seed = 2022, na.rm = FALSE,
                  centre = TRUE, mean = FALSE) {
  if(na.rm == TRUE){
    df <- df[!is.na(df[[index]]),]
  }
  if(class(df[[index]]) == "factor") {
    n <- length(levels(droplevels(df[[index]])))
    group <- levels(droplevels(df[[index]]))
  } else {
    n <- length(unique(df[[index]][!is.na(df[[index]])]))
    group <- sort(unique(df[[index]][!is.na(df[[index]])]))
  }
  tmp <- vector("list",n)
  for (i in seq_len(n)) {
    x <- df[[PC1]][df[[index]] == group[i]]
    y <- df[[PC2]][df[[index]] == group[i]]
    base <- data.frame(x = x, y = y)
    wt <- df[[weight]][df[[index]] == group[i]]
    tmp[[i]]<-gt.test(base = base, wt = wt, centre = centre, mean = mean)
    names(tmp)[i] <- group[i]
  }
  return(tmp)
}

# Function for extracting weighted mean of components/cos2 by category
m_point <- function(df, PC1 = "PC3", PC2 = "PC1", 
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

## Structured Data Analysis: c_alphan, class
# 1. Find mean points for c_alphan
c_wt_mean <- m_point(ISSP5, index = "c_alphan")
# Calculate mean cos^2
c_wt_cos2 <- m_point(ISSP5, "cos2_3", "cos2_1", index = "c_alphan")
c_wt_cos2$xy <- c_wt_cos2$x + c_wt_cos2$y
c_wt_cos2 %>%
  arrange(desc(x))
c_wt_cos2 %>%
  arrange(desc(y))
c_wt_cos2 %>%
  arrange(desc(xy))

# 1.a Find mean point for China (ISSP4)
cn_wt_mean <- m_point(ISSP4, index = "c_alphan")
cn_wt_cos2 <- m_point(ISSP4, "cos2_3", "cos2_1", index = "c_alphan")

# 2. Find mean points for class (v61)
class_wt_mean <- m_point(ISSP5, index = "v61", na.rm = TRUE)
# Calculate mean cos^2
class_wt_cos2 <- m_point(ISSP5, "cos2_3", "cos2_1", 
                         index = "v61", na.rm = TRUE)
class_wt_cos2$xy <- class_wt_cos2$x + class_wt_cos2$y
class_wt_cos2 %>%
  arrange(desc(x))
class_wt_cos2 %>%
  arrange(desc(y))
class_wt_cos2 %>%
  arrange(desc(xy))

# 3. Find mean points for educations (DEGREE_c)
edu_wt_mean <- m_point(ISSP5, index = "DEGREE_c", na.rm = TRUE)
# Calculate mean cos^2
edu_wt_cos2 <- m_point(ISSP5, "cos2_3", "cos2_1", 
                       index = "DEGREE_c", na.rm = TRUE)
edu_wt_cos2 %>%
  arrange(desc(x))
edu_wt_cos2 %>%
  arrange(desc(y))

# 4. Find mean points for ISCO08 classification 
ISSP5$ISCO08_class <- substr(as.character(ISSP5$ISCO08_code),1,1)
isco_wt_mean <- m_point(ISSP5, index = "ISCO08_class", na.rm = TRUE)
# Calculate mean cos^2
isco_wt_cos2 <- m_point(ISSP5, "cos2_3", "cos2_1",
                        index = "ISCO08_class", na.rm = TRUE)
isco_wt_cos2 %>%
  arrange(desc(x))
isco_wt_cos2 %>%
  arrange(desc(y))

# 5. Find mean points for gender
sex_wt_mean <- m_point(ISSP5, index = "SEX", na.rm = TRUE)
sex_wt_cos2 <- m_point(ISSP5, "cos2_3", "cos2_1", index = "SEX", na.rm = TRUE)

# 6. Find mean points for political party
pol_wt_mean <- m_point(ISSP5, index = "PARTY_LR", na.rm = TRUE)
pol_wt_cos2 <- m_point(ISSP5, "cos2_3", "cos2_1", index = "PARTY_LR", na.rm = TRUE)

## Supplementary variables: Soc is, Soc ought
# 1. Find mean points for Type of Society is (v48)
ISSP5$v48_t <- ISSP5$v48
levels(ISSP5$v48_t)[c(1:5)] <- c("A", "B", "C", "D", "E")
socis_wt_mean <- m_point(ISSP5, index = "v48_t", na.rm = TRUE)
# Calculate mean cos^2
socis_wt_cos2 <- m_point(ISSP5, "cos2_3", "cos2_1",
                         index = "v48_t", na.rm = TRUE)
socis_wt_cos2 %>%
  arrange(desc(x))
socis_wt_cos2 %>%
  arrange(desc(y))

## Geometric Typicality Test (each c_alphan compared with origin)
# a. Using mean of principal components (i.e. zero)
c_wt_gt <- by_gt(ISSP5, index = "c_alphan")
saveRDS(c_wt_gt, "data/c_wt_gt_1_3")

# a.ii. China (ISSP4)
cn_wt_gt <- by_gt(ISSP4, index = "c_alphan")
saveRDS(cn_wt_gt, "data/cn_wt_gt_1_3")

# b. Using barycentre of weighted cloud as mean
c_wt_cf_gt <- by_gt(ISSP5, index = "c_alphan",
                     centre = FALSE,
                     mean = c(weighted.mean(ISSP5$PC3, w = ISSP5$WEIGHT_c),
                              weighted.mean(ISSP5$PC1, w = ISSP5$WEIGHT_c)))
saveRDS(c_wt_cf_gt, "data/c_wt_cf_gt_1_3")

## Visualise
## Plot compatibility ellipses (not NZ, RU or US)
# drop NZ, RU and US
c_wt_mean_sig <- c_wt_mean[-c(10,12,15),]
c_wt_gt_sig <- c_wt_gt[-c(10,12,15)]
# data for country ellipses
el_c_wt <- vector("list", length(c_wt_gt_sig))
for(i in 1:length(el_c_wt)) {
  names(el_c_wt) <- names(c_wt_gt_sig)
  el_c_wt[[i]] <- el.data(c_wt_gt_sig[[i]][[4]])
}
el_c_wt$values <- data.frame(
  "e" = do.call("rbind", lapply(el_c_wt, "[[", 3)),
  "theta" = do.call("rbind", lapply(el_c_wt, "[[", 4))
)
# China
el_cn_wt <- el.data(cn_wt_gt[[1]][[4]])
# m_distance circle
radius = 0.4
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) 
# mean of weighted cloud
wt_x <- weighted.mean(ISSP5$PC3,ISSP5$WEIGHT_c)
wt_y <- weighted.mean(ISSP5$PC1,ISSP5$WEIGHT_c)

# function for scaling cos2 to plot alpha shading for PC2
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
c_wt_cos2$xs <- round((range01(c_wt_cos2$x)*10))+1
# function for scaling text colour
colfunc <- colorRampPalette(c("grey70", "black"))
cols <- colfunc(11)
c_wt_cos2$col <- cols[c_wt_cos2$xs]
sig_col <- c_wt_cos2$col[-c(10,12,15)]

# Plot
png(filename = "plots/c_wt_geom_1_3.png",
    width = 1200, height = 1200)
par(pty = "s", cex = 2)
plot(x = NULL,
     y = NULL,
     main = NULL,
     xlab = paste0("PC3, ", round(mod$Vaccounted[[2,3]]*100, 1), "% of variance\n (Capital Structure/Symbolic System)"),
     ylab = paste0("PC1, ", round(mod$Vaccounted[[2,1]]*100, 1), "% of variance\n (Capital Volume)"),
     xlim = c(-.8,.8),
     ylim = c(-.8,1.15),
     asp = 1)
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y,
      lty = 1, col = "grey80")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
text(c(1.05,1.05,-1.05,-1.05,0,0),c(0.03,-0.03,0.03,-0.03,1.15,-0.8),
     pos = c(2,2,4,4,3,1),
     label = c("INHERITANCE", "CULTURAL CAPITAL+",
               "SOLIDARITY", "SOCIAL CAPITAL+",
               "CAPITAL VOLUME+", "CAPITAL VOLUME-"),
     cex = 0.8)
# Add text
text(c_wt_mean$x+c(-0.1,0,0,0,0,-0.05,0,0,0.05,rep(0,7)), 
     c_wt_mean$y+c(0,0,0,0,0,0.08,0,0,-0.08,rep(0,7)),
     label = c_wt_mean$group,
     col = c_wt_cos2$col)
# China
text(cn_wt_mean$x, cn_wt_mean$y-0.07,
     label = "CN (2009)", col = cols[1])
# Compatibility region for mean
for(i in 1:length(c_wt_gt_sig)) {
  ellipse(c(c_wt_mean_sig$x[i],c_wt_mean_sig$y[i]),
          c_wt_gt_sig[[i]][[4]],
          c_wt_gt_sig[[i]][[5]],
          center.pch = FALSE,
          col = sig_col[i],
          lwd = 1)
}
# China
ellipse(c(cn_wt_mean$x, cn_wt_mean$y),
        cn_wt_gt[[1]][[4]],
        cn_wt_gt[[1]][[5]],
        center.pch = FALSE,
        col = cols[1])
# Vertices
for(i in 1:length(c_wt_gt_sig)) {
  vertices.k(c(c_wt_mean_sig$x[i], 
               c_wt_mean_sig$y[i]),
             c_wt_gt_sig[[i]][[4]],
             c_wt_gt_sig[[i]][[5]],
             col = sig_col[i])
}
# China
vertices.k(c(cn_wt_mean$x, cn_wt_mean$y),
           cn_wt_gt[[1]][[4]],
           cn_wt_gt[[1]][[5]],
           col = "grey60")
# Mean point of weighted cloud
points(wt_x,wt_y,pch=19,col="grey90",cex=2)
# Mean points class
points(class_wt_mean$x,
       class_wt_mean$y,
       pch = 15, col = "black")
lines(class_wt_mean$x[2:5], class_wt_mean$y[2:5], 
      lty = 4)
text(class_wt_mean$x, class_wt_mean$y, pos = c(3,2,4,1,4,3),
     label = c("lower","working","lower middle","middle","upper middle","upper"))
# Mean points of Society Is/Ought
points(socis_wt_mean$x,socis_wt_mean$y,
       pch = 17, col = "grey50")
text(socis_wt_mean$x,socis_wt_mean$y, pos = c(4,4,4,4,2),
     label = socis_wt_mean$group, col = "grey50")
lines(socis_wt_mean$x[1:4],socis_wt_mean$y[1:4],
      lty = 4)
legend("topleft", legend = c("Social Class","Type of Soc. Is",
                              "Weighted Barycenter"),
       col = c("black","grey50","grey90"), pch = c(15,17,19), cex = 1)
dev.off()

## Geometric Typicality Test (each class compared with origin)
class_wt_gt <- by_gt(ISSP5, index = "v61", na.rm = TRUE)
saveRDS(class_wt_gt, "data/class_wt_gt_1_3")

## Visualise
## Plot compatibility ellipse for lower class
# m_distance circle
radius = 0.4
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) 

## Plot
png(filename = "results/class_wt_geom_1_3.png",
    width = 1200, height = 1200)
par(pty = "s", cex = 2)
plot(x = NULL,
     y = NULL,
     main = "Class Mean Points and Compatibility Ellipses, \n PC1 and PC3",
     xlab = paste0("PC3, ", round(mod$Vaccounted[[2,3]]*100, 1), "% of variance\n (Symbolic System)"),
     ylab = paste0("PC1, ", round(mod$Vaccounted[[2,1]]*100, 1), "% of variance\n (Capital Volume)"),
     xlim = c(-.8,.8),
     ylim = c(-.8,.8),
     asp = 1)
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y,
      lty = 1, col = "grey80")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# Add other classes as points
points(class_wt_mean$x,
       class_wt_mean$y,
       pch = 13, col = "black")
text(class_wt_mean$x, class_wt_mean$y+0.03,
     label = c("lower","working","lower middle","middle","upper middle","upper"))
# Mean point for weighted cloud
points(wt_x,wt_y,pch=19,col="grey95",cex=2)
# Mean points of Society Is/Ought
points(socis_wt_mean$x,socis_wt_mean$y,
       pch = 17, col = "grey50")
text(socis_wt_mean$x,socis_wt_mean$y+0.03,
     label = socis_wt_mean$group, col = "grey50")
legend("topright", legend = c("Type of Soc. Is",
                              "Weighted Barycenter"),
       col = c("grey50","grey90"), pch = c(17,19), cex = 1)
dev.off()

## Geometric Typicality Test (each ISCO08 class)
isco_wt_gt <- by_gt(ISSP5, index = "ISCO08_class", na.rm = TRUE)
saveRDS(isco_wt_gt, "data/isco_wt_gt_1_3")

# m_distance circle
radius = 0.4
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200)

## Plot
png(filename = "results/isco_wt_geom_1_3.png",
    width = 1200, height = 1200)
par(pty = "s", cex = 2)
plot(x = NULL,
     y = NULL,
     main = "ISCO08 Class Mean Points \n PC1 and PC3",
     xlab = paste0("PC3, ", round(mod$Vaccounted[[2,3]]*100, 1), "% of variance\n (Symbolic System)"),
     ylab = paste0("PC1, ", round(mod$Vaccounted[[2,1]]*100, 1), "% of variance\n (Capital Volume)"),
     xlim = c(-.8,.8),
     ylim = c(-.8,.8),
     asp = 1)
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y,
      lty = 1, col = "grey80")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# Add other classes as points
points(isco_wt_mean$x,isco_wt_mean$y,
       pch = 13, col = "black")
text(isco_wt_mean$x, isco_wt_mean$y-0.03,
     label = c("Manager","Professional","Technician","Clerical",
               "Services","Agri","Craft","Industrial","Elementary"))
# Mean point of weighted cloud
points(wt_x,wt_y,pch=19,col="grey95",cex=2)
# Mean points of Society Is/Ought
points(socis_wt_mean$x,socis_wt_mean$y,
       pch = 17, col = "grey50")
text(socis_wt_mean$x,socis_wt_mean$y+0.03,
     label = socis_wt_mean$group, col = "grey50")
legend("topright", legend = c("Type of Soc. Is",
                              "Weighted Barycenter"),
       col = c("grey50","black"), pch = c(17,19), cex = 1)
dev.off()

## Geometric Typicality Test (each edu compared with origin)
edu_wt_gt <- by_gt(ISSP5, index = "DEGREE_c", na.rm = TRUE)
saveRDS(edu_wt_gt, "data/edu_wt_gt_1_3")

# m_distance circle
radius = 0.4
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) 

## Plot
png(filename = "results/edu_wt_geom_1_3.png",
    width = 1200, height = 1200)
par(pty = "s", cex = 2)
plot(x = NULL,
     y = NULL,
     main = "Education Mean Points, PC1 and PC3",
     xlab = paste0("PC3, ", round(mod$Vaccounted[[2,3]]*100, 1), "% of variance\n (Symbolic System)"),
     ylab = paste0("PC1, ", round(mod$Vaccounted[[2,1]]*100, 1), "% of variance\n (Capital Volume)"),
     xlim = c(-.8,.8),
     ylim = c(-.8,.8),
     asp = 1)
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y,
      lty = 1, col = "grey80")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# Mean point of weighted cloud
points(wt_x,wt_y,pch=19,col="grey95",cex=2)
# Add classes as points
points(edu_wt_mean$x,edu_wt_mean$y,
       pch = 13, col = "black")
text(edu_wt_mean$x, edu_wt_mean$y+0.03,
     label = c("no formal","primary","lower secondary","upper secondary",
               "post secondary","lower tertiary","upper tertiary"))
# Mean points of Society Is/Ought
points(socis_wt_mean$x,socis_wt_mean$y,
       pch = 17, col = "grey50")
text(socis_wt_mean$x,socis_wt_mean$y+0.03,
     label = socis_wt_mean$group, col = "grey50")
legend("topright", legend = c("Type of Soc. Is",
                              "Weighted Barycenter"),
       col = c("grey50","grey95"), pch = c(17,19), cex = 1)
dev.off()

## Geometric Typicality Test (each sex compared with origin)
sex_wt_gt <- by_gt(ISSP5, index = "SEX", na.rm = TRUE)
saveRDS(sex_wt_gt, "data/sex_wt_gt_1_3")

# m_distance circle
radius = 0.4
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) 

## Plot
png(filename = "results/sex_wt_geom_1_3.png",
    width = 1200, height = 1200)
par(pty = "s", cex = 2)
plot(x = NULL,
     y = NULL,
     main = "Gender Mean Points in PCs 1 and 3",
     xlab = paste0("PC3, ", round(mod$Vaccounted[[2,3]]*100, 1), "% of variance\n (Capital Structure)"),
     ylab = paste0("PC1, ", round(mod$Vaccounted[[2,1]]*100, 1), "% of variance\n (Capital Volume)"),
     xlim = c(-.8,.8),
     ylim = c(-.8,.8),
     asp = 1)
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y,
      lty = 1, col = "grey80")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# Mean point of weighted cloud
points(wt_x,wt_y,pch=19,col="grey95",cex=2)
# Add classes as points
points(sex_wt_mean$x,sex_wt_mean$y,
       pch = 13, col = "black")
text(sex_wt_mean$x, sex_wt_mean$y+0.03,
     label = sex_wt_mean$group)
# Mean points of Society Is/Ought
points(socis_wt_mean$x,socis_wt_mean$y,
       pch = 17, col = "grey50")
text(socis_wt_mean$x,socis_wt_mean$y+0.03,
     label = socis_wt_mean$group, col = "grey50")
legend("topright", legend = c("Type of Soc. Is",
                              "Weighted barycenter"),
       col = c("grey50","grey95"), pch = c(17,19), cex = 1)
dev.off()

## Geometric Typicality Test (each pol compared with origin)
pol_wt_gt <- by_gt(ISSP5, index = "PARTY_LR", na.rm = TRUE)
saveRDS(pol_wt_gt, "data/pol_wt_gt_1_3")

# m_distance circle
radius = 0.4
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) 

## Plot
png(filename = "results/pol_wt_geom_1_3.png",
    width = 1200, height = 1200)
par(pty = "s", cex = 2)
plot(x = NULL,
     y = NULL,
     main = "Politics Mean Points in PCs 1 and 3",
     xlab = paste0("PC3, ", round(mod$Vaccounted[[2,3]]*100, 1), "% of variance\n (Capital Structure)"),
     ylab = paste0("PC1, ", round(mod$Vaccounted[[2,1]]*100, 1), "% of variance\n (Capital Volume)"),
     xlim = c(-.8,.8),
     ylim = c(-.8,.8),
     asp = 1)
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y,
      lty = 1, col = "grey80")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# Mean point of weighted cloud
points(wt_x,wt_y,pch=19,col="grey95",cex=2)
# Add classes as points
points(pol_wt_mean$x,pol_wt_mean$y,
       pch = 13, col = "black")
text(pol_wt_mean$x, pol_wt_mean$y+0.03,
     label = pol_wt_mean$group)
# Mean points of Society Is/Ought
points(socis_wt_mean$x,socis_wt_mean$y,
       pch = 17, col = "grey50")
text(socis_wt_mean$x,socis_wt_mean$y+0.03,
     label = socis_wt_mean$group, col = "grey50")
legend("topright", legend = c("Type of Soc. Is",
                              "Weighted barycenter"),
       col = c("grey50","grey95"), pch = c(17,19), cex = 1)
dev.off()