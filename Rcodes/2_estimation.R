rm(list = ls())

# --------- library and functions -------------
source("Rcodes/0_library.R")
source("Rcodes/0a_helper.R")

# --------- load data ------------------------
load("data/meta-data-generation-seed1.RData")

# -------- estimation, corhort only ----------
test = mem.mono(y = cohort.dat$y, 
                              z = cohort.dat$z, 
                              va = cbind(cohort.dat$v1, cohort.dat$v2),
                              vb = cbind(cohort.dat$v1, cohort.dat$v2),
                              alpha.start = c(0,0),
                              beta.start = c(0,0))

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
ps = exp(predict(ps.model)) / (exp(predict(ps.model)) + 1)
ps.mat = cbind(1-ps, ps)

test3 = max.likelihood.casecrt(y = casecrt.dat$y, 
                               z = casecrt.dat$z, 
                               va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                               vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                               ps.mat = ps.mat,
                               alpha.start = c(0,0),
                               beta.start = c(0,0))
test4 = max.likelihood.casecrt(y = casecrt.dat$y, 
                               z = casecrt.dat$z, 
                               va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                               vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                               ps.mat = ps.mat,
                               alpha.start = test$point.est[c(1,2)],
                               beta.start = test$point.est[c(3,4)])
# -------- estimation, both ----------
ps.model = glm(c(cohort.dat$y, casecrt.dat$z) ~ c(cohort.dat$v2,casecrt.dat$z), family = binomial(link='logit'))

ps = exp(predict(ps.model)) / (exp(predict(ps.model)) + 1)
ps.mat = cbind(1-ps, ps)[-c(1:nrow(cohort.dat)), ]

test.all = max.likelihood.all(case.control = list(y = casecrt.dat$y, 
                                                  z = casecrt.dat$z, 
                                                  va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                                  vb = cbind(casecrt.dat$v1, casecrt.dat$v2)),
                              cohort = list(y = cohort.dat$y, 
                                            z = cohort.dat$z, 
                                            va = cbind(cohort.dat$v1, cohort.dat$v2),
                                            vb = cbind(cohort.dat$v1, cohort.dat$v2),
                                            alpha.start = c(0,0),
                                            beta.start = c(0,0)),
                              ps.mat = ps.mat,
                              alpha.start = c(0,0),
                              beta.start = c(0,0))
test.all$coefficients

test$coefficients
test1$coefficients

test3$coefficients
test4$coefficients

test$coefficients
