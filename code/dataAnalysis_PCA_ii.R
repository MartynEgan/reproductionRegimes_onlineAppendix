######################
# Data Analysis: PCA #
######################
# Analyse model #
#################

# Load packages
library("dplyr")
library("psych")

# Load data
ISSP5 <- readRDS("data/ISSP5_GDA")
ISSP4 <- readRDS("data/ISSP4_GDA")
mod <- readRDS("data/pca_mod")
scores <- readRDS("data/pca_scores")

# Analysis of components
# Check sums, column space
sum(mod$loadings[]^2) # sum of squares of all loadings equals number of vars
near(colSums(mod$loadings^2), mod$values) # sum of squares of loadings per component equals eigenvalue
rowSums(mod$loadings^2) # sum of squares per variable equals 1
sum(mod$values) # eigenvalues sum to number of variables

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
colSums(contribs) # columns sum to 1
rowSums(contribs) # rows sum to 1

# Biplot of components
# Note: requires edited version of biplot.psych function, 
# biplot.edit(), found in plotone_edit.R
source("code/plotoneEdit.R")

# A) First and second principal components
# 1. create plot data object (x = PC2, y = PC1)
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
text_x$loadings[1,2] <- text_x$loadings[1,2] - 0.03
text_x$loadings[1,1] <- text_x$loadings[1,1] + 0.03
text_x$loadings[5,2] <- text_x$loadings[5,2] - 0.06
text_x$loadings[5,1] <- text_x$loadings[5,1] + 0.03

# 3. some code for plotting unit circle https://stackoverflow.com/questions/22265704/drawing-circle-in-r
radius = 1
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) 

# 4. plot 
png("plots/biplot_PCA.png", res = 300, height = 1600, width = 1600)
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

# B) First and third principal components
# 1. create plot data object
x <- list()
x$scores <- scores$scores[,c(3,1)] 
x$scores[,2] <- x$scores[,2] *-1 # change sign of PC1 scores to reflect direction of underlying Likert scale
x$loadings <- mod$loadings[,c(3,1)]
x$loadings[,1] <- x$loadings[,1] *-1 # change sign of PC3 loadings to match direction of PC2 in previous plot
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
text_x$loadings[2,2] <- text_x$loadings[2,2]-0.03

# 3. some code for plotting unit circle https://stackoverflow.com/questions/22265704/drawing-circle-in-r
radius = 1
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) 

# 4. plot 
png("plots/biplot_PCA_PC3.png", res = 300, height = 1600, width = 1600)
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

saveRDS(ISSP5, "data/ISSP5_GDA")

for (i in 1:length(coords)) {
  ISSP4[,paste0("cos2_",i)] <- ISSP4[,coords[i]]^2 / (rowSums(ISSP4[,coords]^2))
}

saveRDS(ISSP4, "data/ISSP4_GDA")

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
