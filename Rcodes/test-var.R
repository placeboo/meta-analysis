rm(list = ls())
library(tidyverse)
library(brm)

n1 = 1000
n2 = 500

source("Rcodes/0a_helper.R")
source("Rcodes/0b_variance.R")
source("Rcodes/0c_simulation.R")

gamma.true = c(0,-1)
beta.true = c(-7, 0.5)
eta.true = c(0.2,-1)

# --------- load data ------------------------
cohort.dat = simCohort(n = n1, gamma.true = gamma.true, beta.true = beta.true, eta.true = eta.true)
casectrl.dat = simCasecrt(n = n2, pr = 1/3, gamma.true = gamma.true, beta.true = beta.true, eta.true = eta.true)

test1 = brm(y = cohort.dat$y, 
            x = cohort.dat$z, 
            va = cbind(cohort.dat$v1, cohort.dat$v2),
            vb = cbind(cohort.dat$v1, cohort.dat$v2),
            param = "RR",
            alpha.start = c(0,0),
            beta.start = c(0,0))

ps.model = glm(z ~ v2, data = cohort.dat, family = binomial(link='logit'))
ps = predict(ps.model, type = "response", newdata = casectrl.dat)
ps.mat = cbind(1-ps,ps)
test.all = max.likelihood.all(case.control = list(y = casectrl.dat$y, 
                                                  z = casectrl.dat$z, 
                                                  va = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                                  vb = cbind(casectrl.dat$v1, casectrl.dat$v2)),
                              cohort = list(y = cohort.dat$y, 
                                            z = cohort.dat$z, 
                                            va = cbind(cohort.dat$v1, cohort.dat$v2),
                                            vb = cbind(cohort.dat$v1, cohort.dat$v2)),
                              alpha.start = test1$point.est[c(1,2)],
                              beta.start = test1$point.est[c(3,4)],
                              eta.start = ps.model$coefficients)

# variance function -------------------
var.mat = var.meta(gamma = test.all$point.est[c(1,2)],
         beta = test.all$point.est[c(3,4)],
         eta = test.all$point.est[c(5,6)],
         case.control =list(y = casectrl.dat$y, 
                            z = casectrl.dat$z, 
                            va = cbind(casectrl.dat$v1, casectrl.dat$v2),
                            vb = cbind(casectrl.dat$v1, casectrl.dat$v2),
                            vc = cbind(casectrl.dat$v1, casectrl.dat$v2)),
         cohort = list(y = cohort.dat$y, 
                                z = cohort.dat$z, 
                                va = cbind(cohort.dat$v1, cohort.dat$v2),
                                vb = cbind(cohort.dat$v1, cohort.dat$v2),
                                vc = cbind(cohort.dat$v1, cohort.dat$v2)))


sqrt(diag(var.mat)) / test.all$se.est


var.casectrl.mat = var.casectrl(gamma = test.all$point.est[c(1,2)],
                            beta = test.all$point.est[c(3,4)],
                            eta = test.all$point.est[c(5,6)],
                            case.control =list(y = casectrl.dat$y, 
                                               z = casectrl.dat$z, 
                                               va = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                               vb = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                               vc = cbind(casectrl.dat$v1, casectrl.dat$v2)),
                            cohort =  list(z = cohort.dat$z,
                                           va = cbind(cohort.dat$v1, cohort.dat$v2), 
                                           vb = cbind(cohort.dat$v1, cohort.dat$v2), 
                                           vc = cbind(cohort.dat$v1, cohort.dat$v2)))
 
sqrt(diag(var.casectrl.mat)) 


# case-control only variance comparison ------
test.casectrl = max.likelihood.casectrl(y = casectrl.dat$y, 
                                        z = casectrl.dat$z, 
                                        va = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                        vb = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                        ps.mat = ps.mat,
                                        alpha.start = c(0,0),
                                        beta.start = c(0,0))

test.casectrl2 = max.likelihood.casectrl2(y = casectrl.dat$y, 
                                          z = casectrl.dat$z, 
                                          va = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                          vb = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                          alpha.start = c(0,0),
                                          beta.start = c(0,0),
                                          eta.start =  ps.model$coefficients)

test.casectrl3 = max.likelihood.casectrl3(case.control = list(y = casectrl.dat$y, 
                                                              z = casectrl.dat$z, 
                                                              va = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                                              vb = cbind(casectrl.dat$v1, casectrl.dat$v2)),
                                          cohort = list(z = cohort.dat$z, 
                                                        va = cbind(cohort.dat$v1, cohort.dat$v2),
                                                        vb = cbind(cohort.dat$v1, cohort.dat$v2)),
                                          alpha.start = test1$point.est[c(1,2)],
                                          beta.start = test1$point.est[c(3,4)],
                                          eta = ps.model$coefficients)


test.casectrl$se.est
test.casectrl2$se.est
test.casectrl3$se.est

test.casectrl3$se.est / test.casectrl2$se.est

fct.se = sqrt(diag(var.casectrl(gamma = test.casectrl$point.est[c(1,2)],
                                beta = test.casectrl$point.est[c(3,4)],
                                eta = ps.model$coefficients,
                                case.control =list(y = casectrl.dat$y, 
                                                   z = casectrl.dat$z, 
                                                   va = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                                   vb = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                                   vc = cbind(casectrl.dat$v1, casectrl.dat$v2)),
                                cohort = list(z = cohort.dat$z, 
                                              va = cbind(cohort.dat$v1, cohort.dat$v2),
                                              vb = cbind(cohort.dat$v1, cohort.dat$v2),
                                              vc = cbind(cohort.dat$v1, cohort.dat$v2)))))

fct.se3 = sqrt(diag(var.casectrl(gamma = test.casectrl3$point.est[c(1,2)],
                                beta = test.casectrl3$point.est[c(3,4)],
                                eta = test.casectrl3$point.est[c(5,6)],
                                case.control =list(y = casectrl.dat$y, 
                                                   z = casectrl.dat$z, 
                                                   va = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                                   vb = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                                   vc = cbind(casectrl.dat$v1, casectrl.dat$v2)),
                                cohort = list(z = cohort.dat$z, 
                                              va = cbind(cohort.dat$v1, cohort.dat$v2),
                                              vb = cbind(cohort.dat$v1, cohort.dat$v2),
                                              vc = cbind(cohort.dat$v1, cohort.dat$v2)))))

fct.se3/test.casectrl3$se.est

