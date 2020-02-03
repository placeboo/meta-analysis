rm(list = ls())

# --------- library and functions -------------
source("Rcodes/0_library.R")

source("Rcodes/0a_helper.R")

# --------- load data ------------------------
n1 = 1000
n2 = 500
seed = 2
load( paste("data/crt", n1, "_casecrt", n2, "_seed", seed, ".RData", sep = ""))

test1 = brm(y = cohort.dat$y, 
            x = cohort.dat$z, 
            va = cbind(cohort.dat$v1, cohort.dat$v2),
            vb = cbind(cohort.dat$v1, cohort.dat$v2),
            param = "RR",
            alpha.start = c(0,0),
            beta.start = c(0,0))



# -------- estimation, case control only ----------
# propensity score model
ps.model = glm(casecrt.dat$z ~ casecrt.dat$v2, family = binomial(link='logit'))


test2 = max.likelihood.casecrt(y = casecrt.dat$y, 
                               z = casecrt.dat$z, 
                               va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                               vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                               alpha.start = c(0,0),
                               beta.start = c(0,0),
                               eta.start = c(0,0))

test3 = max.likelihood.casecrt(y = casecrt.dat$y, 
                               z = casecrt.dat$z, 
                               va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                               vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                               alpha.start = c(0,0),
                               beta.start = c(0,0),
                               eta.start = ps.model$coefficients)
# -------- estimation, both ----------
ps.model1 = glm(c(cohort.dat$z, casecrt.dat$z) ~ c(cohort.dat$v2, casecrt.dat$v2), family = binomial(link='logit'))


test4 = max.likelihood.casecrt(y = casecrt.dat$y, 
                               z = casecrt.dat$z, 
                               va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                               vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                               alpha.start = test1$point.est[c(1,2)],
                               beta.start = test1$point.est[c(3,4)],
                               eta.start = ps.model1$coefficients)

test.all = max.likelihood.all(case.control = list(y = casecrt.dat$y, 
                                                  z = casecrt.dat$z, 
                                                  va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                                  vb = cbind(casecrt.dat$v1, casecrt.dat$v2)),
                              cohort = list(y = cohort.dat$y, 
                                            z = cohort.dat$z, 
                                            va = cbind(cohort.dat$v1, cohort.dat$v2),
                                            vb = cbind(cohort.dat$v1, cohort.dat$v2)),
                              alpha.start = c(0,0),
                              beta.start = c(0,0),
                              eta.start = ps.model1$coefficients)
                              



test1$coefficients
test2$coefficients
test3$coefficients
test4$coefficients
test.all$coefficients

