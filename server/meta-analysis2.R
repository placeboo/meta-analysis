#!/usr/local/bin/Rscript
rm(list = ls()) 

args = commandArgs(TRUE)

n1 = as.numeric(args[[1]])
seed = as.numeric(args[[2]])

library(brm)
library(dplyr)
source("Rcodes/0a_helper.R")
source("Rcodes/0c_simulation.R")


# ------ study setting ------------
st = "02082020"

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)
n2.vec = c(200, 500, 1000, 2000, 4000)


do.one.simulation = function(n1, n2.vec, pr, gamma.true, beta.true, eta.true) {
     
     cohort.dat = simCohort(n1, gamma.true, beta.true, eta.true)
     
     test1 = brm(y = cohort.dat$y, 
                 x = cohort.dat$z, 
                 va = cbind(cohort.dat$v1, cohort.dat$v2),
                 vb = cbind(cohort.dat$v1, cohort.dat$v2),
                 param = "RR",
                 alpha.start = c(0,0),
                 beta.start = c(0,0))
     
     
     
     # -------- estimation, case control only ----------
     # propensity score model
     ps.model = glm(z ~ v2, data = cohort.dat, family = binomial(link='logit'))
     
     opt = c(test1$coefficients[,c(1,2)])
     
     for (n2 in n2.vec) {
          
          casecrt.dat = simCasecrt(n2, pr, gamma.true, beta.true, eta.true)     
          ps = predict(ps.model, type = "response", newdata = casecrt.dat)
          ps.mat = cbind(1-ps, ps)
          
          test2 = max.likelihood.casecrt(y = casecrt.dat$y, 
                                         z = casecrt.dat$z, 
                                         ps.mat = ps.mat,
                                         va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                         vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                         alpha.start = test1$point.est[c(1,2)],
                                         beta.start = test1$point.est[c(3,4)])
          
          # -------- estimation, both ----------
          test.all = max.likelihood.all(case.control = list(y = casecrt.dat$y, 
                                                            z = casecrt.dat$z, 
                                                            va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                                            vb = cbind(casecrt.dat$v1, casecrt.dat$v2)),
                                        cohort = list(y = cohort.dat$y, 
                                                      z = cohort.dat$z, 
                                                      va = cbind(cohort.dat$v1, cohort.dat$v2),
                                                      vb = cbind(cohort.dat$v1, cohort.dat$v2)),
                                        alpha.start = test1$point.est[c(1,2)],
                                        beta.start = test1$point.est[c(3,4)],
                                        eta.start = ps.model$coefficients)
          
          opt = c(opt, c(test2$coefficients[,c(1,2)]), c(test.all$coefficients[, c(1,2)]))
     }
      
     return(opt)
}

set.seed(seed)
N_sim = 20
est.mat = replicate(N_sim, do.one.simulation(n1, n2.vec, pr = 1/3,  gamma.true, beta.true, eta.true))

index = unique(unlist(apply(est.mat, 1, function(x) which(is.na(x)))))

while (length(index) > 0){
     est.mat = est.mat[,-index]
     est.tmp = replicate(length(index), do.one.simulation(n1, n2, n21, gamma.true, beta.true, eta.true))
     est.mat = cbind(est.mat, est.tmp)
     
     index = unique(unlist(apply(est.mat, 1, function(x) which(is.na(x)))))
     
}


save(est.mat, file = paste("data6/crt", n1, "_casecrt", n2, "_seed", seed, "_", st, ".RData", sep = ""))



