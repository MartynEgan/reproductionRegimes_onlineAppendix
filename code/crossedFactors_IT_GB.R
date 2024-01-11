##################################################
# Data Analysis: Crossed Factors (IT and GB-GBN) #
##################################################

#######
# PCA #
#######

# Load packages
library("dplyr")
library("psych")

# Load data
ISSP5 <- readRDS("data/ISSP5_x_drop")

# Filter IT and GB-GBN
ISSP5 <- ISSP5 %>%
  filter(c_alphan %in% c("IT","GB-GBN"))

# Active variables
a_vars_5 <- c("num_v1", # wealthy family
              "num_v2", # well-educated parents
              "num_v3", # own education
              "num_v4", # hard work
              "num_v5", # knowing right people
              "num_v6", # political connections
              "num_v7") # bribes

# Run PCA with weights and no rotation
mod <- pca(ISSP5[,a_vars_5],
           nfactors = 7,
           cor = "poly",
           rotate = "none",
           weight = ISSP5$WEIGHT_c,
           scores = FALSE)

# Derive scores 
scores <- factor.scores(ISSP5[,a_vars_5],
                        f = mod)

# Bind scores
ISSP5 <- cbind(ISSP5, scores$scores)

# Analyse components
round(mod$Vaccounted[1:3,],3) # eigenvalues and variance explained
round(mod$loadings[],3) # correlation of variables to each component
round(mod$loadings[]^2,3) # Cos^2 of each variable to each axis
contribs <- matrix(rep(NA, 7*7), nrow = 7, ncol = 7,
                   dimnames = list(row.names(mod$loadings),
                                   colnames(mod$loadings)))
for (i in 1:ncol(mod$loadings)) {
  contribs[,i] <- mod$loadings[,i]^2 / mod$values[i] 
}
round(contribs, 3) # contribution of variables to inertia of axis

# Biplot of components
# Note: requires edited version of biplot.psych function, 
# biplot.edit(), found in plotoneEdit.R
source("code/plotoneEdit.R")

# A) First and second principal components
# 1. create plot data object
x <- list()
x$scores <- scores$scores[,c(2,1)]
x$scores[,2] <- x$scores[,2] *-1 # change sign of PC1 scores to reflect direction of underlying Likert scale
x$scores[,1] <- x$scores[,1] *-1 # change sign of PC2 scores to reflect true direction of component
x$loadings <- mod$loadings[,c(2,1)]
row.names(x$loadings) <- c("EC_fam",
                           "CC_fam",
                           "CC",
                           "hard_work",
                           "SC",
                           "SC_pol",
                           "bribes")
class(x) <- c("psych","fa")

# 2. edit text positions
text_x <- x
text_x$loadings[1,2] <- text_x$loadings[1,2]-0.07
text_x$loadings[4,2] <- text_x$loadings[4,2]+0.03
text_x$loadings[5,1] <- text_x$loadings[5,1]+0.04
text_x$loadings[5,2] <- text_x$loadings[5,2]-0.03

# 3. some code for plotting unit circle https://stackoverflow.com/questions/22265704/drawing-circle-in-r
radius = 1
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) 

# 4. plot 
png("plots/biplot_PCA_IT_GB.png", res = 300, height = 1600, width = 1600)
biplot.edit(x,
            xlim.s = c(-3,3),
            ylim.s = c(-3,3),
            main = NULL,
            cex = c(0.7,0.7),
            col = c("grey", "black"),
            smoother = TRUE,
)
text(x = text_x$loadings[,1], 
     y = text_x$loadings[,2], 
     labels = row.names(x$loadings), 
     cex = 0.8, 
     col = "black")
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y)
dev.off()

# Analysis of scores
# Check sums, row space
coords <- ISSP5 %>%
  select(starts_with("PC")) %>%
  colnames()

round(sapply(ISSP5[,coords], mean),2) # mean of scores per component = 0
round(sapply(ISSP5[,coords], sd),2) # std. dev. of scores per component = 0
round(cor(ISSP5[,coords], use = "pairwise"),2) # PCs are uncorrelated in row space

# Cos^2 of individuals on each axis
for (i in 1:length(coords)) {
  ISSP5[,paste0("cos2_",i)] <- ISSP5[,coords[i]]^2 / (rowSums(ISSP5[,coords]^2))
}

saveRDS(ISSP5, "data/IT_GB_GDA")

cos2 <- ISSP5 %>%
  select(starts_with("cos2")) %>%
  colnames()

rowSums(ISSP5[,cos2]) # cos2 sums to 1 across all rows

# Find mean cos^2 for each country across first 3 axes
country_cos2 <- aggregate(ISSP5$cos2_1, by = list(ISSP5$c_alphan), FUN = mean)
country_cos2$x2 <- aggregate(ISSP5$cos2_2, by = list(ISSP5$c_alphan), FUN = mean)$x
country_cos2$x3 <- aggregate(ISSP5$cos2_3, by = list(ISSP5$c_alphan), FUN = mean)$x
country_cos2$x <- round(country_cos2$x, 3)
country_cos2$x2 <- round(country_cos2$x2, 3)
country_cos2$x3 <- round(country_cos2$x3, 3)
country_cos2 <- country_cos2 %>%
  arrange(desc(x))
country_cos2

#######
# GDA #
#######

# Change sign of PC1 and PC2 scores to reflect direction of Likert scale
ISSP5$PC1 <- ISSP5$PC1 *-1
ISSP5$PC2 <- ISSP5$PC2 *-1

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
  theta <- atan(E[[2]][[2,1]]/E[[2]][[1,1]]) * (180/pi)
  out <- list()
  out$eig <- E
  out$D <- D
  out$e <- e
  out$theta <- theta
  return(out)
}

# Function for performing geometric typicality test
source("code/gt_test.R")
by_gt <- function(df, PC1 = "PC2", PC2 = "PC1", index, 
                  weight = "WEIGHT_c", notable_D = 0.4, 
                  alpha = 0.05, max_number = 1e+06, 
                  n_dir = 200, seed = 2022, na.rm = FALSE) {
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
    tmp[[i]]<-gt.test(base = base, wt = wt)
    names(tmp)[i] <- group[i]
  }
  return(tmp)
}

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

## Structured Data Analysis: c_alphan, class
# 1. Find mean points for c_alphan
c_wt_mean <- m_point(ISSP5, index = "c_alphan")
# Calculate mean cos^2
c_wt_cos2 <- m_point(ISSP5, "cos2_2", "cos2_1", index = "c_alphan")
c_wt_cos2 %>%
  arrange(desc(x))
c_wt_cos2 %>%
  arrange(desc(y))

# 2. Find mean points for class (v61)
class_wt_mean <- m_point(ISSP5, index = "v61", na.rm = TRUE)
# Calculate mean cos^2
class_wt_cos2 <- m_point(ISSP5, "cos2_2", "cos2_1", 
                         index = "v61", na.rm = TRUE)
class_wt_cos2 %>%
  arrange(desc(x))
class_wt_cos2 %>%
  arrange(desc(y))

# 3. Find mean points for class by country
IT_class_wt_mean <- m_point(ISSP5[ISSP5$c_alphan == "IT",], 
                            index = "v61", na.rm = TRUE)
IT_class_wt_cos2 <- m_point(ISSP5[ISSP5$c_alphan == "IT",], 
                            "cos2_2", "cos2_1", 
                            index = "v61", na.rm = TRUE)
IT_class_wt_cos2 %>%
  arrange(desc(x))
IT_class_wt_cos2 %>%
  arrange(desc(y))

GB_class_wt_mean <- m_point(ISSP5[ISSP5$c_alphan == "GB-GBN",], 
                            index = "v61", na.rm = TRUE)
GB_class_wt_cos2 <- m_point(ISSP5[ISSP5$c_alphan == "GB-GBN",], 
                            "cos2_2", "cos2_1", 
                            index = "v61", na.rm = TRUE)
GB_class_wt_cos2 %>%
  arrange(desc(x))
GB_class_wt_cos2 %>%
  arrange(desc(y))

## Supplementary variable: Soc is
# 1. Find mean points for Type of Society is (v48)
ISSP5$v48_t <- ISSP5$v48
levels(ISSP5$v48_t)[c(1:5)] <- c("A", "B", "C", "D", "E")
socis_wt_mean <- m_point(ISSP5, index = "v48_t", na.rm = TRUE)
# Calculate mean cos^2
socis_wt_cos2 <- m_point(ISSP5, "cos2_2", "cos2_1",
                         index = "v48_t", na.rm = TRUE)
socis_wt_cos2 %>%
  arrange(desc(x))
socis_wt_cos2 %>%
  arrange(desc(y))

## Geometric Typicality Test (each c_alphan compared with origin)
c_wt_gt <- by_gt(ISSP5, index = "c_alphan")
saveRDS(c_wt_gt, "data/c_wt_gt_IT_GB")

## Visualise
## Plot compatibility ellipses 
# data for country ellipses
el_c_wt <- vector("list", length(c_wt_gt))
for(i in 1:length(el_c_wt)) {
  names(el_c_wt) <- names(c_wt_gt)
  el_c_wt[[i]] <- el.data(c_wt_gt[[i]][[4]])
}
el_c_wt$values <- data.frame(
  "e" = do.call("rbind", lapply(el_c_wt, "[[", 3)),
  "theta" = do.call("rbind", lapply(el_c_wt, "[[", 4))
)

# m_distance circle
radius = 0.4
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) 
# mean of weighted cloud
wt_x <- weighted.mean(ISSP5$PC2,ISSP5$WEIGHT_c)
wt_y <- weighted.mean(ISSP5$PC1,ISSP5$WEIGHT_c)

# Plot
png(filename = "plots/cf_wt_IT_GB.png",
    width = 1200, height = 1200)
par(pty = "s", cex = 2)
plot(x = NULL,
     y = NULL,
     main = NULL,
     xlab = paste0("PC2, ", round(mod$Vaccounted[[2,2]]*100, 1), "% of variance\n (Capital Structure/Symbolic System)"),
     ylab = paste0("PC1, ", round(mod$Vaccounted[[2,1]]*100, 1), "% of variance\n (Capital Volume)"),
     xlim = c(-0.8,0.8),
     ylim = c(-0.8,0.8),
     asp = 1)
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y,
      lty = 1, col = "grey80")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# add class points
points(class_wt_mean[class_wt_mean$group != c("Lower class", "Upper class"),"x"],
       class_wt_mean[class_wt_mean$group != c("Lower class", "Upper class"),"y"],
       pch = 15, col = "black")
points(IT_class_wt_mean[IT_class_wt_mean$group != c("Lower class", "Upper class"),"x"],
       IT_class_wt_mean[IT_class_wt_mean$group != c("Lower class", "Upper class"),"y"],
       pch = 0, col = "black")
text(IT_class_wt_mean$x[2:5], IT_class_wt_mean$y[2:5],
     pos = c(2,1,4,4),
     label = c("working","lower\nmiddle","middle","upper middle"))
points(GB_class_wt_mean[GB_class_wt_mean$group != c("Lower class", "Upper class"),"x"],
       GB_class_wt_mean[GB_class_wt_mean$group != c("Lower class", "Upper class"),"y"],
       pch = 7, col = "black")
text(GB_class_wt_mean$x[2:5], GB_class_wt_mean$y[2:5],
     pos = c(1,1,4,1),
     label = c("working","lower\nmiddle","middle","upper\nmiddle"))
# add lines, within country
lines(IT_class_wt_mean$x[2:5], IT_class_wt_mean$y[2:5], lty = 1)
lines(GB_class_wt_mean$x[2:5], GB_class_wt_mean$y[2:5], lty = 1)
# add lines, within class
for(i in 2:5){
lines(c(IT_class_wt_mean$x[i],
        class_wt_mean$x[i],
        GB_class_wt_mean$x[i]),
      c(IT_class_wt_mean$y[i],
        class_wt_mean$y[i],
        GB_class_wt_mean$y[i]), lty = 2)
}
legend("topright", legend = c("IT","GB","Mean"),
       col = c("black","black","black"), pch = c(0,7,15), cex = 1)
dev.off()

saveRDS(ISSP5, "data/IT_GB_GDA")

#############
# Inference #
#############

# 1.Prepare data
ISSP5 <- readRDS("data/IT_GB_GDA")
ISSP5 <- ISSP5[!is.na(ISSP5$v61),]
ISSP5 <- ISSP5[!(ISSP5$v61 %in% c("Lower class", "Upper class")),]
ISSP5$v61 <- droplevels(ISSP5$v61)
Y.IL <- as.matrix(ISSP5[, c("PC2","PC1")]) # principal coordinates
n <- dim(Y.IL)[1] # number of cases
L <- dim(Y.IL)[2] # number of principal axes
Lambda.L <- colSums(Y.IL^2)/n # variances of axes
A <- as.factor(ISSP5[, "c_alphan"]) # A factor (c_alphan)
B <- ISSP5[, "v61"] # B factor (Social Class)
C <- factor(paste(as.vector(A),as.vector(B), sep = "")) # crossing AxB
n.C <- as.vector(table(C)) # frequencies of groups
names(n.C) <- levels(C); cardC <- length(n.C) # number of groups
n.C_G <- as.vector(table(A))
names(n.C_G) <- levels(A); cardA <- length(n.C_G)
n.Cl <- as.vector(table(B))
names(n.Cl) <- levels(B); cardB <- length(n.Cl)
#C_int.C <- matrix(c(1,-1,-1,1), nrow = 1) # interaction of contrast and
C_int.C <- matrix(c(1,-1,-1,1,1,-1,-1,1),nrow = 1)
w_int <- 1/(n * sum(n.C^-1)) # inverse of its squared norm

# 2.Descriptive analysis
# coordinates of mean points
ISSP5$C <- C
ISSP5$A <- A
ISSP5$B <- B
Y_G.CL <- as.matrix(m_point(ISSP5, index = "C")[,-1])
#Y_G.CL <- as.matrix(aggregate(Y.IL, by = list(C), FUN = mean)[,-1])
C_G.CL <- as.matrix(m_point(ISSP5, index = "A")[,-1])
#C_G.CL <- as.matrix(aggregate(Y.IL, by = list(A), FUN = mean)[,-1])
Cl.CL <- as.matrix(m_point(ISSP5, index = "B")[,-1])
#Cl.CL <- as.matrix(aggregate(Y.IL, by = list(B), FUN = mean)[,-1])
cat(" Weights and means of mean points,\n")
print(cbind(n.C,Y_G.CL))
print(cbind(n.C_G,C_G.CL))
print(cbind(n.Cl, Cl.CL))
# variances of sources of variation
V_AxB.L <- n.C %*% Y_G.CL^2/n
V_A.L <- n.C_G %*% C_G.CL^2/n
V_B.L <- n.Cl %*% Cl.CL^2/n
V_int.L <- w_int * (C_int.C %*% Y_G.CL)^2 # try orthogonal contrast matrix with 5x6 matrix
V_add.L <- V_AxB.L - V_int.L
row.names(V_AxB.L) <- "AxB"
row.names(V_A.L) <- "A"
row.names(V_B.L) <- "B"
row.names(V_int.L) <- "A.B"
row.names(V_add.L) <- "A+B"
# printing of results
cat("\n Decomposition of variances\n")
tmp1 <- rbind(V_AxB.L, V_add.L, V_int.L, V_A.L, V_B.L)
tmp2 <- t(t(c(sum(V_AxB.L), sum(V_add.L), sum(V_int.L), sum(V_A.L), sum(V_B.L))))
colnames(tmp2) <- "Plane 1_2"
tmp3 <- round((cbind(tmp1,tmp2)), digits=5)
A_in_B <- tmp3[1,] - tmp3[5,]
B_in_A <- tmp3[1,] - tmp3[4,]
tmp3 <- rbind(tmp3, A_in_B, B_in_A)
print(tmp3)
cf_var <- tmp3
eta_2 <- sum(V_AxB.L)/sum(Lambda.L)

# # step 3. Inductive analysis
# # observed M-variances in plane 1-2
VM_AxB <- sum(V_AxB.L %*% diag(1/Lambda.L))
VM_int <- sum(V_int.L %*% diag(1/Lambda.L))
VM_add <- VM_AxB - VM_int
VM_A <- sum(V_A.L %*% diag(1/Lambda.L))
VM_B <- sum(V_B.L %*% diag(1/Lambda.L))
# distributions of test statistics
Z.IL <- Y.IL %*% diag(1/sqrt(Lambda.L)) # standardised principal coordinates
cardJ <- 1000000 # number of nestings
VM.AxB.J <- VM_int.J <- VM_add.J <- VM_A.J <- VM_B.J <- matrix(0, cardJ)
Z_Cj.CL <- matrix(0, nrow=cardC, ncol = L)
for(j in 1:cardJ) {
  AA <- sample(C)
  for (c in (1:cardC)) {    # coordinates of mean points for J
    one.C <- rep(0L, cardC); one.C[c]=1
    Z_Cj.CL[c, ] <- t(one.C[AA]) %*% Z.IL/n.C[c]
  }
  # M-variances of sources
  VM.AxB.J[j] <- sum(n.C %*% Z_Cj.CL^2)/n     # AxB
  VM_int.J[j] <- w_int * sum((C_int.C %*% Z_Cj.CL)^2) # A.B
}
VM_add.J <- VM.AxB.J - VM_int.J                 # A+B

n_AxB <- sum(VM.AxB.J >= VM_AxB*(1-1e-12)) # testing AxB

n_int <- sum(VM_int.J >= VM_int*(1-1e-12)) # testing A.B and A+B
n_add <- sum(VM_add.J >= VM_add*(1-1e-12))
cat("\n","p-value\n for A+B ", n_add,"/", cardJ, " = ",
    round(n_add/cardJ, digits = 3), "\n for A.B ", n_int,"/", cardJ, " = ",
    round(n_int/cardJ, digits = 3), sep = "")

p_AxB <- round(n_AxB/cardJ, digits = 3)
p_int <- round(n_int/cardJ, digits = 3)
p_add <- round(n_add/cardJ, digits = 3)

# Main Effects
VM_A.J <- matrix(0, cardJ)
Z_Cj.C_G <- matrix(0, nrow=cardA, ncol = L)
for(j in 1:cardJ) {
  AA <- sample(A)
  for (c in (1:cardA)) {    # coordinates of mean points for J
    one.A <- rep(0L, cardA); one.A[c]=1
    Z_Cj.C_G[c, ] <- t(one.A[AA]) %*% Z.IL/n.C_G[c]
  }
  # M-variances of sources
  VM_A.J[j] <- sum(n.C_G %*% Z_Cj.C_G^2)/n     # A
}

n_A <- sum(VM_A.J >= VM_A*(1-1e-12)) # testing A
cat("\n","p-value\n for A ", n_A,"/", cardJ, " = ",
    round(n_A/cardJ, digits = 3))
p_A <- round(n_A/cardJ, digits = 3)

VM_B.J <- matrix(0, cardJ)
Z_Cj.Cl <- matrix(0, nrow=cardB, ncol = L)
for(j in 1:cardJ) {
  AA <- sample(B)
  for (c in (1:cardB)) {    # coordinates of mean points for J
    one.B <- rep(0L, cardB); one.B[c]=1
    Z_Cj.Cl[c, ] <- t(one.B[AA]) %*% Z.IL/n.Cl[c]
  }
  # M-variances of sources
  VM_B.J[j] <- sum(n.Cl %*% Z_Cj.Cl^2)/n     # B
}

n_B <- sum(VM_B.J >= VM_B*(1-1e-12)) # testing B
cat("\n","p-value\n for B ", n_B,"/", cardJ, " = ",
    round(n_B/cardJ, digits = 3))
p_B <- round(n_B/cardJ, digits = 3)

p_vals <- matrix(c(p_A,p_B,p_AxB,p_int,p_add),
                 dimnames = list(c("p_A","p_B","p_AxB","p_int","p_add"),c("p_value")),
                 ncol=1)
cf_IT_GB <- list()
cf_IT_GB$cf_var <- cf_var
cf_IT_GB$eta_2 <- eta_2
cf_IT_GB$p_vals <- p_vals
cf_IT_GB$Lambda <- Lambda.L
saveRDS(cf_IT_GB, "data/cf_IT_GB")