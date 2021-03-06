#!/usr/local/bin/Rscript
rm(list = ls()) 

args = commandArgs(TRUE)

n1 = as.numeric(args[[1]])
seed = as.numeric(args[[2]])

library(brm)
library(dplyr)

source("Rcodes/0b_helper_1Dim.R") 


# ------ study setting ------------
st = "02082020"

gamma.true = c(0.5)
beta.true = c(-7)
eta.true = c(0.2)
n2.vec = c(200, 500, 1000, 2000, 4000)

# ------- data simulation function------

simCohort = function(n, gamma.true, beta.true, eta.true) {
     
     va = as.matrix(rep(1, n))
     pz = exp(va %*% eta.true) / (1 + exp(va %*% eta.true))
     z = rbinom(n, 1, pz)
     
     P.mat = t(mapply(getProbScalarRR, va %*% gamma.true, va %*% beta.true))
     
     y = rep(0, n)
     
     for(i in 1:2){
          y[z == i-1] = rbinom(length(which(z==i-1)), 1, P.mat[z==i-1,i])
     }
     
     data = as.data.frame(cbind(va, z, y))
     names(data) = c("va", "z", "y")
     
     return(data)
     
}

simCasecrt = function(n, pr, gamma.true, beta.true, eta.true) {
     n1 = round(n * pr) 
     n0 = n - n1
     
     data = simCohort(n, gamma.true, beta.true, eta.true)
     
     control = data %>%
          filter(y == 0)
     
     case = data %>%
          filter(y == 1)
     
     if ((nrow(control) < n0) & (nrow(case) > n1)) {
          
          case = case[1: n1, ]
          while(nrow(control) < n0) {
               tmp.dat = simCohort(n, gamma.true, beta.true, eta.true)
               control = tmp.dat %>%
                    filter(y == 0) %>%
                    rbind(control)
          }
          
          control = control[1:n0, ]
          
     } else if ((nrow(control) > n0) & (nrow(case) < n1)) {
          control = control[1: n0, ]
          
          while(nrow(case) < n1) {
               tmp.dat =simCohort(n, gamma.true, beta.true, eta.true)
               
               case = tmp.dat %>%
                    filter(y == 1) %>%
                    rbind(case)
          }
     }
     
     case = case[1:n1,]
     
     return(rbind(case, control))
     
}
do.one.simulation = function(n1, n2.vec, pr, gamma.true, beta.true, eta.true) {
     
     #------------- meta dataset generation ------------
     
     cohort.dat = simCohort(n1, gamma.true, beta.true, eta.true)
     
     test1 = mem.mono(y = cohort.dat$y,
                      z = cohort.dat$z,
                      va = as.matrix(cohort.dat$va),
                      vb = as.matrix(cohort.dat$va),
                      alpha.start = 0,
                      beta.start = 0)
     
     
     opt = c(test1$coefficients[,c(1,2)])
     ps = mean(cohort.dat$z)
     
     # -------- estimation, both ----------
     for (n2 in n2.vec) {
          
          casecrt.dat = simCasecrt(n2, pr, gamma.true, beta.true, eta.true)   
          
          test.all = max.likelihood.all(case.control = list(y = casecrt.dat$y, 
                                                            z = casecrt.dat$z, 
                                                            va = as.matrix(casecrt.dat$va),
                                                            vb = as.matrix(casecrt.dat$va)),
                                        cohort = list(y = cohort.dat$y, 
                                                      z = cohort.dat$z, 
                                                      va = as.matrix(cohort.dat$va),
                                                      vb = as.matrix(cohort.dat$va)),
                                        alpha.start = test1$point.est[1],
                                        beta.start = test1$point.est[2],
                                        eta.start = ps)
          
          opt = c(opt,
                  c(test.all$coefficients[,c(1,2)]))     
          }

     return(opt)
}

N_sim = 50
set.seed(seed)
est.mat = replicate(N_sim, do.one.simulation(n1, n2.vec, pr = 1/3, gamma.true, beta.true, eta.true))

index = unique(unlist(apply(est.mat, 1, function(x) which(is.na(x)))))

while (length(index) > 0){
     est.mat = est.mat[,-index]
     est.tmp = replicate(length(index), do.one.simulation(n1, n2, n21, gamma.true, beta.true, eta.true))
     est.mat = cbind(est.mat, est.tmp)
     
     index = unique(unlist(apply(est.mat, 1, function(x) which(is.na(x)))))
     
}

save(est.mat, file = paste("data-sanity/crt", n1, "_seed", seed, "_", st, ".RData", sep = ""))



