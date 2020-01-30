#!/usr/local/bin/Rscript
rm(list = ls())

args = commandArgs(TRUE)

n = as.numeric(args[[1]])
seed = as.numeric(args[[2]])

library(brm)

set.seed(seed)
# ------ study setting ------------
st = "01292020"

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)


# ------- data simulation function------

simData = function(n, gamma.true, beta.true, eta.true) {
     
     v1 = rep(1, n)
     v2 = runif(n, -2, 2)
     va = cbind(v1, v2)
     
     pz = exp(va %*% eta.true) / (1 + exp(va %*% eta.true))
     z = rbinom(n, 1, pz)
     
     P.mat = t(mapply(getProbScalarRR, va %*% gamma.true, va %*% beta.true))
     
     y = rep(0, n)
     
     for(i in 1:2){
          y[z == i-1] = rbinom(length(which(z==i-1)), 1, P.mat[z==i-1,i])
     }
     
     data = as.data.frame(cbind(va, z, y))
     names(data) = c("v1", "v2", "z", "y")
     
     return(data)
     
}

do.one.simulation = function(n, gamma.true, beta.true, eta.true) {
     #------------- meta dataset generation ------------
     cohort.dat = simData(n, gamma.true, beta.true, eta.true)
     print(sum(cohort.dat$y)/n)
     #save(cohort.dat, casecrt.dat, file = paste("data/crt", n1, "_casecrt", n2, "_seed", seed, "_",st, ".RData", sep = ""))
     
     # -------- estimation, corhort only ----------
     test1 = brm(y = cohort.dat$y, 
                 x = cohort.dat$z, 
                 va = cbind(cohort.dat$v1, cohort.dat$v2),
                 vb = cbind(cohort.dat$v1, cohort.dat$v2),
                 param = "RR",
                 alpha.start = c(0,0),
                 beta.start = c(0,0))
     

     opt = c(test1$coefficients[,c(1,2)])
     
     return(opt)
}


N_sim = 2
est.mat = replicate(N_sim, do.one.simulation(n, gamma.true, beta.true, eta.true))

save(est.mat, file = paste("crt/crt", n, "_seed", seed, "_", st, ".RData", sep = ""))



