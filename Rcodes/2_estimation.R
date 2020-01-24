rm(list = ls())

# --------- library and functions -------------
source("Rcodes/0_library.R")
source("Rcodes/0a_helper.R")

# --------- load data ------------------------
load("data/meta-data-generation-seed1.RData")

# -------- estimation, corhort only ----------
test = max.likelihood.v2(y = cohort.dat$y, 
                              z = as.factor(cohort.dat$z), 
                              va = cbind(cohort.dat$v1, cohort.dat$v2),
                              vb = cbind(cohort.dat$v1, cohort.dat$v2),
                              alpha.start = c(0,0),
                              beta.start = c(0,0))


library(MEM)
test2 = mem.mono(y = cohort.dat$y, 
         z = cohort.dat$z, 
         va = cbind(cohort.dat$v1, cohort.dat$v2),
         vb = cbind(cohort.dat$v1, cohort.dat$v2),
         alpha.start = c(0,0),
         beta.start = c(0,0))
# -------- estimation, both ----------
