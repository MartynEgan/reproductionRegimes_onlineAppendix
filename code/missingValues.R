###################################
# Investigation of missing values #
###################################

# Load packages
library("dplyr")

# Load data
ISSP4_x <- readRDS("data/ISSP4_x")
ISSP5_x <- readRDS("data/ISSP5_x")

# Active variables
a_vars_5 <- c("v1", # wealthy family
              "v2", # well-educated parents
              "v3", # own education
              "v4", # hard work
              "v5", # knowing right people
              "v6", # political connections
              "v7") # bribes

# Check and missing values in active variables
ISSP4_x %>%
  select(all_of(a_vars_5)) %>%
  filter(if_any(all_of(a_vars_5), ~ is.na(.))) %>%
  summarise(across(all_of(a_vars_5), ~sum(is.na(.) == TRUE)))

ISSP5_x %>%
  select(all_of(a_vars_5)) %>%
  filter(if_any(all_of(a_vars_5), ~ is.na(.))) %>%
  summarise(across(all_of(a_vars_5), ~sum(is.na(.) == TRUE)))

# Note: it appears that the category "can't choose" has been coded NA in
# the dataset. This is a significant number for v7 (bribery), and most 
# likely represents an important category of respondent; unfortunately
# because of this coding error the category cannot be explored in analysis.

# Compare proportions of NAs across countries for v7
v7_table <- left_join(as.data.frame(table(ISSP5_x[,"c_alphan"])),
                      as.data.frame(table(ISSP5_x[is.na(ISSP5_x$v7) == TRUE,"c_alphan"])),
                      by = "c_alphan",
                      suffix = c(".all", ".NA"))
rownames(v7_table) <- v7_table[,1]
v7_table <- v7_table[,2:3]
v7_mean_NA <- sum(v7_table$Freq.NA)/sum(v7_table$Freq.all)
prop_test_v7 <- prop.test(v7_table[,2], v7_table[,1])
prop_test_v7_mean <- list(NULL)
for (i in 1:nrow(v7_table)) {
  prop_test_v7_mean[[i]] <- prop.test(v7_table[i,2],
                                      v7_table[i,1],
                                      p = v7_mean_NA,
                                      alternative = "greater")
}
v7_table$p_val <- round(do.call("rbind", lapply(prop_test_v7_mean, "[[", 3)),2)

# Germany, Finland, France, Great Britain, Japan, Norway, Russia and Sweden
# all have proportions of missing values for v7 significantly greater than
# the average.

# Compare proportions of NAs across all active variables
v_NA_table <- ISSP5_x %>%
  select(all_of(a_vars_5)) %>%
  filter(if_any(all_of(a_vars_5), ~ is.na(.))) %>%
  summarise(across(all_of(a_vars_5), ~sum(is.na(.) == TRUE)))

v_NA_table <- t(v_NA_table)
v_NA_table <- as.data.frame(v_NA_table)
v_NA_table[,2] <- rep(nrow(ISSP5_x), nrow(v_NA_table))

prop.test(v_NA_table[,1],
          v_NA_table[,2])

mean_v_NA <- mean(v_NA_table[,1]) / nrow(ISSP5_x)

prop_test_v_NA_mean <- list(NULL)
for (i in 1:nrow(v_NA_table)) {
  prop_test_v_NA_mean[[i]] <- prop.test(v_NA_table[i,1],
                                        v_NA_table[i,2],
                                        p = mean_v_NA,
                                        alternative = "greater")
}
v_NA_table$p_val <- round(do.call("rbind", lapply(prop_test_v_NA_mean, "[[", 3)),2)

# Both v6 and v7 feature missing value counts which are significantly
# different to the mean of missing values across active variables.

# All NAs in active variables are dropped.
ISSP4_x_drop <- ISSP4_x[complete.cases(ISSP4_x[,a_vars_5]) == TRUE,]
ISSP5_x_drop <- ISSP5_x[complete.cases(ISSP5_x[,a_vars_5]) == TRUE,]

# Compare distributions on key variables.
# by country
round(prop.table(table(ISSP5_x$c_alphan))*100, 2)
round(prop.table(table(ISSP5_x_drop$c_alphan))*100, 2)
# by social class
round(prop.table(table(ISSP5_x$v61))*100, 2)
round(prop.table(table(ISSP5_x_drop$v61))*100, 2)