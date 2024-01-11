######################
# Data Analysis: PCA #
######################
# Run model #
#############

# Load packages
library("dplyr")
library("psych")

# Load data
ISSP5 <- readRDS("data/ISSP5_x_drop")
ISSP4 <- readRDS("data/ISSP4_x_drop")

# Active variables
a_vars_5 <- c("num_v1", # wealthy family
              "num_v2", # well-educated parents
              "num_v3", # own education
              "num_v4", # hard work
              "num_v5", # knowing right people
              "num_v6", # political connections
              "num_v7") # bribes

# Compare Pearson and Polychoric correlation matrices
prs <- ISSP5 %>%
  select(all_of(a_vars_5)) %>%
  as.matrix() %>%
  cor(method = "pearson")

pco <- ISSP5 %>%
  select(all_of(a_vars_5)) %>%
  as.matrix() %>%
  polychoric() 

cor.plot(prs) 
cor.plot(pco$rho)

# correlations are improved using polychoric correlation

# Perform PCA on Polychoric correlation matrix
# a) Test for optimal components
ISSP5 %>%
  select(all_of(a_vars_5)) %>%
  as.matrix() %>%
  fa.parallel(cor = "poly", fm = "minres", fa = "pc", 
              sim = FALSE)

# b) Run PCA with weights and no rotation
mod <- pca(ISSP5[,a_vars_5],
           nfactors = 7,
           cor = "poly",
           rotate = "none",
           weight = ISSP5$WEIGHT_c,
           scores = FALSE)

# c) Derive scores 
# ISSP5
scores <- factor.scores(ISSP5[,a_vars_5],
                        f = mod)
# ISSP4 (passive individuals)
p_scores <- factor.scores(ISSP4[,a_vars_5],
                          f = mod)

# Bind scores
ISSP5 <- cbind(ISSP5, scores$scores)
ISSP4 <- cbind(ISSP4, p_scores$scores)

# Save results
saveRDS(ISSP5, "data/ISSP5_GDA")
saveRDS(ISSP4, "data/ISSP4_GDA")
saveRDS(mod, "data/pca_mod")
saveRDS(scores, "data/pca_scores")
