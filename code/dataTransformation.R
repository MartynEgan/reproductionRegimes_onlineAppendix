##############################
# Data Transformation Script #
##############################

library("dplyr")

# Load data
ISSP5 <- readRDS("data/ZA7600v3.RDS")
ISSP4 <- readRDS("data/ZA5400.RDS")

# Active variables, structuring factors and passive variables
a_vars_4 <- c("V6", # wealthy family
              "V7", # well-educated parents
              "V8", # own education
              "V10",# hard work
              "V11",# knowing right people
              "V12",# political connections
              "V13")# bribes

a_vars_5 <- c("v1", # wealthy family
              "v2", # well-educated parents
              "v3", # own education
              "v4", # hard work
              "v5", # knowing right people
              "v6", # political connections
              "v7") # bribes

s_fact_4 <- c("C_ALPHAN", # Country abbreviation
              "V66", # Social class (self-reported)
              "AGE", # Age (recode)
              "SEX", # Sex
              "DEGREE", # Highest education level - different levels to ISSP4
              "WRKST", # Main status (employment) - different levels to MAINSTAT
              "V69", # Father's Occupation
              "V70", # Mother's Occupation
              "V72", # Own Occupation
              "ISCO88", # Occupation ISCO 1988
              "PARTY_LR", # Party affiliation left-right 
              "URBRURAL", # Type of community (self-reported)
              "WEIGHT" # Sampling weights
              )
              
s_fact_5 <- c("c_alphan", # Country abbreviation
              "v61", # Social class (self-reported)
              "AGE", # Age (recode)
              "SEX", # Sex of respondent
              "DEGREE", # Highest completed degree (derived) - different levels to ISSP4
              "MAINSTAT", # Main status (employment) - different levels to WRKST
              "v67", # Father's Occupation
              "v68", # Mother's Occupation
              "v69", # Own occupation
              "v70", # Partner's Occupation
              "ISCO08", # Occupation ISCO 2008
              "PARTY_LR", # Party voted for last election
              "URBRURAL", # Place of living: urban - rural
              "WEIGHT" # Sampling weights
              )

p_vars_4 <- c("V54", # Type of country is
              "V55", # Type of country prefer
              "V32" # Differences in income too large - different to v50
              )
  
p_vars_5 <- c("v48", # Type of country is
              "v49", # Type of country prefer
              "v50", # How fair is income distribution - different to V32
              "v53" # How difficult to make ends meet
              )

# Country selection
ctry_4 <- c("AU", # Australia
            "GB-GBN", # Great Britain
            "US", # United States
            "NZ", # New Zealand
            "JP", # Japan
            "TW", # Taiwan
            "DE", # Germany
            "FR", # France
            "IT", # Italy
            "SE", # Sweden
            "DK", # Denmark
            "FI", # Finland
            "NO", # Norway
            "RU", # Russia
            "VE", # Venezuela
            "PH", # Philippines
            "CN" # China (missing in ISSP5)
            )
  
ctry_5 <- c("AU", # Australia
            "GB-GBN", # Great Britain
            "US", # United States
            "NZ", # New Zealand
            "JP", # Japan
            "TW", # Taiwan
            "DE", # Germany
            "FR", # France
            "IT", # Italy
            "SE", # Sweden
            "DK", # Denmark
            "FI", # Finland
            "NO", # Norway
            "RU", # Russia
            "VE", # Venezuela
            "PH" # Philippines
            )

# Filter on rows and columns
ISSP4_x <- filter(ISSP4, C_ALPHAN %in% ctry_4)
ISSP5_x <- filter(ISSP5, c_alphan %in% ctry_5)

ISSP4_x <- ISSP4_x %>%
  select(all_of(c(a_vars_4,
                  s_fact_4,
                  p_vars_4)))
ISSP5_x <- ISSP5_x %>%
  select(all_of(c(a_vars_5,
                  s_fact_5,
                  p_vars_5)))

# Rename ISSP4 variables to match ISSP5
ISSP4_x <- ISSP4_x %>%
  rename_with(~ a_vars_5, all_of(a_vars_4)) %>%
  rename_with(~ s_fact_5[!s_fact_5 %in% "v70"], all_of(s_fact_4)) %>%
  rename_with(~ p_vars_5[-4], all_of(p_vars_4))

# Relevel ISSP4 active variables to match ISSP5
# Check for unused levels and drop
sapply(ISSP4_x[,1:7], table)
sapply(ISSP5_x[,1:7], table)
ISSP4_x[,1:7] <- droplevels(ISSP4_x[,1:7])
ISSP5_x[,1:7] <- droplevels(ISSP5_x[,1:7])

# Check and drop levels for other variables
sapply(ISSP4_x[,c(9,11,12,13,14,15,16,18,19,21,22,23)], table)
sapply(ISSP5_x[,c(9,11,12,13,14,15,16,17,19,20,22,23,24)], table)
ISSP4_x[,c(9,11,12,13,14,15,16,18,19,21,22,23)] <- droplevels(ISSP4_x[,c(9,11,12,13,14,15,16,18,19,21,22,23)])
ISSP5_x[,c(9,11,12,13,14,15,16,17,19,20,22,23,24)] <- droplevels(ISSP5_x[,c(9,11,12,13,14,15,16,17,19,20,22,23,24)])

# Convert AGE to categorical variable
ISSP4_x$AGE_c <- ISSP4_x$AGE
levels(ISSP4_x$AGE_c)[levels(ISSP4_x$AGE_c)=="98 years or more"] <- "98"
ISSP4_x$AGE_c <- as.numeric(as.character(ISSP4_x$AGE_c))
ISSP4_x$AGE_c <- cut(ISSP4_x$AGE_c, 
                     breaks = c(-Inf, 35, 50, 65, Inf),
                     labels = c("15-34", "35-49", "50-64", "65+"),
                     right = FALSE)

ISSP5_x$AGE_c <- ISSP5_x$AGE
levels(ISSP5_x$AGE_c)[levels(ISSP5_x$AGE_c)=="15 years"] <- "15"
levels(ISSP5_x$AGE_c)[levels(ISSP5_x$AGE_c)=="89 years; US: 89 years or more"] <- "89"
levels(ISSP5_x$AGE_c)[levels(ISSP5_x$AGE_c)=="105 years"] <- "105"
ISSP5_x$AGE_c <- as.numeric(as.character(ISSP5_x$AGE_c))
ISSP5_x$AGE_c <- cut(ISSP5_x$AGE_c, 
                     breaks = c(-Inf, 35, 50, 65, Inf),
                     labels = c("15-34", "35-49", "50-64", "65+"),
                     right = FALSE)

# Tidy levels for DEGREE
ISSP5_x$DEGREE_c <- ISSP5_x$DEGREE
levels(ISSP5_x$DEGREE_c) <- c("No formal education",
                              "Primary School",
                              "Lower Secondary",
                              "Upper Secondary",
                              "Post secondary",
                              "Lower tertiary",
                              "Upper tertiary")

# Save intermediate data objects (see investigation
# of missing values script)

saveRDS(ISSP4_x, "data/ISSP4_x")
saveRDS(ISSP5_x, "data/ISSP5_x")

# All NAs in active variables are dropped.
ISSP4_x_drop <- ISSP4_x[complete.cases(ISSP4_x[,a_vars_5]) == TRUE,]
ISSP5_x_drop <- ISSP5_x[complete.cases(ISSP5_x[,a_vars_5]) == TRUE,]

# Weight ISSP5 cases by n of country
country_n <- table(ISSP5_x_drop$c_alphan)
base <- 1000
c_weight <- base / country_n
c_weight <- as.data.frame(c_weight)
names(c_weight) <- c("c_alphan", "weight")

ISSP5_x_drop <- left_join(ISSP5_x_drop, c_weight)
ISSP5_x_drop$WEIGHT_c <- ISSP5_x_drop$WEIGHT * ISSP5_x_drop$weight 
# Align weight to sum to nrow
constant <- nrow(ISSP5_x_drop)/sum(ISSP5_x_drop$WEIGHT_c)
ISSP5_x_drop$WEIGHT_c <- constant*ISSP5_x_drop$WEIGHT_c

# Repeat for ISSP4
country_n <- table(ISSP4_x_drop$c_alphan)
base <- 1000
c_weight <- base / country_n
c_weight <- as.data.frame(c_weight)
names(c_weight) <- c("c_alphan", "weight")

ISSP4_x_drop <- left_join(ISSP4_x_drop, c_weight)
ISSP4_x_drop$WEIGHT_c <- ISSP4_x_drop$WEIGHT * ISSP4_x_drop$weight 
# Align weight to sum to nrow
constant <- nrow(ISSP4_x_drop)/sum(ISSP4_x_drop$WEIGHT_c)
ISSP4_x_drop$WEIGHT_c <- constant*ISSP4_x_drop$WEIGHT_c

# Pre-processing: convert active variables to numeric
ISSP5_x_drop <- ISSP5_x_drop %>%
  mutate(across(.cols = all_of(a_vars_5),
                .fns = as.numeric,
                .names = "num_{.col}"))

ISSP4_x_drop <- ISSP4_x_drop %>%
  mutate(across(.cols = all_of(a_vars_5),
                .fns = as.numeric,
                .names = "num_{.col}"))

# Check cast to numeric for both data sets
head(ISSP4_x_drop[ISSP4_x_drop$v1 == "Essential","num_v1"])
head(ISSP5_x_drop[ISSP5_x_drop$v1 == "Essential","num_v1"])

# Prop table for active vars
round(prop.table(table(ISSP5_x$v1, useNA = "ifany"))*100,1)

# Save data
saveRDS(ISSP4_x_drop, "data/ISSP4_x_drop")
saveRDS(ISSP5_x_drop, "data/ISSP5_x_drop")