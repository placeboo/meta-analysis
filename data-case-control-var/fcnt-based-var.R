#!/usr/local/bin/Rscript
rm(list = ls()) 

args = commandArgs(TRUE)

n = as.numeric(args[[1]])
seed = as.numeric(args[[2]])

library(brm)
library(dplyr)
source("Rcodes/0a_helper.R")
source("Rcodes/0b_variance.R")
source("Rcodes/0c_simulation.R")

# ------ study setting ------------
st = "02242020"

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)

do.one.simulation = function(n, pr, gamma.true, beta.true, eta.true) {
     # cohort have eta
     cohort.dat = simCohort(n = 1000, gamma.true, beta.true, eta.true)
     
     ps.model = glm(cohort.dat$z ~ cohort.dat$v2, family = binomial(link='logit'))
     
     casectrl.dat = simCasecrt(n, pr, gamma.true, beta.true,  eta.true)
     
     test.casectrl2 = max.likelihood.casectrl2(y = casectrl.dat$y, 
                                               z = casectrl.dat$z, 
                                               va = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                               vb = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                               alpha.start = c(0,0),
                                               beta.start = c(0,0), 
                                               eta.start = ps.model$coefficients)
     
     fct.se2 = sqrt(diag(var.casectrl(gamma = test.casectrl2$point.est[c(1,2)],
                                      beta = test.casectrl2$point.est[c(3,4)],
                                      eta = ps.model$coefficients,
                                      case.control =list(y = casectrl.dat$y, 
                                                         z = casectrl.dat$z, 
                                                         va = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                                         vb = cbind(casectrl.dat$v1, casectrl.dat$v2),
                                                         vc = cbind(casectrl.dat$v1, casectrl.dat$v2)))))
     
     
     opt = c(test.casectrl2$coefficients[, c(1,2)], fct.se2)
     
     return(opt)
}

set.seed(seed)
N_sim = 200

est.mat = replicate(N_sim, do.one.simulation(n, pr=1/3, gamma.true, beta.true, eta.true))

index = unique(unlist(apply(est.mat, 1, function(x) which(is.na(x)))))

while (length(index) > 0){
     est.mat = est.mat[,-index]
     est.tmp = replicate(length(index), do.one.simulation(n, pr=1/3, gamma.true, beta.true, eta.true))
     est.mat = cbind(est.mat, est.tmp)
     
     index = unique(unlist(apply(est.mat, 1, function(x) which(is.na(x)))))
     
}


save(est.mat, file = paste("data-case-control-var/casecrt", n, "_seed", seed, "_", st, ".RData", sep = ""))
