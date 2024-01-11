########################################
# File conversion script for ZA7600 v3 #
########################################

# Load packages
library("haven") 
library("dplyr")

# Load data (note: must be acquired from GESIS website)
ZA7600 <- read_spss("data/ZA7600_v3-0-0.sav")

# Transform labels into factors for all other variables
ZA7600$WEIGHT <- as.numeric(ZA7600$WEIGHT)
ZA7600 <- haven::as_factor(ZA7600)

# Save as RDS
saveRDS(ZA7600, "data/ZA7600v3.RDS")

#################################
# File conversion script ZA5400 #
#################################

# Load data (note: must be acquired from GESIS website)
ZA5400 <- read_spss("data/ZA5400_v4-0-0.sav")

# Transform labels into factors for all other variables
ZA5400$WEIGHT <- as.numeric(ZA5400$WEIGHT)
ZA5400 <- haven::as_factor(ZA5400)

# Save as RDS
saveRDS(ZA5400, "data/ZA5400.RDS")