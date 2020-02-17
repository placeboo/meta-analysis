#!/usr/local/bin/Rscript
rm(list = ls())

args = commandArgs(TRUE)

n1 = as.numeric(args[[1]])
n2 = as.numeric(args[[2]])
seed = as.numeric(args[[3]])

library(brm)
library(dplyr)
source("0a_helper.R")

set.seed(seed)
# ------ study setting ------------
st = "01292020"
n21 = round(n2/3)

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)


# ------- data simulation function------

simData = function(n1, n2, n21, gamma.true, beta.true, eta.true) {
     
     # cohort study
        v1 = rep(1, 20000)
        v2 = runif(20000, -2, 2)
        va = cbind(v1, v2)
        
        pz = exp(va %*% eta.true) / (1 + exp(va %*% eta.true))
        z = rbinom(20000, 1, pz)
        
        P.mat = t(mapply(getProbScalarRR, va %*% gamma.true, va %*% beta.true))
        
        y = rep(0, 20000)
        
        for(i in 1:2){
                y[z == i-1] = rbinom(length(which(z==i-1)), 1, P.mat[z==i-1,i])
        }
        
        data = as.data.frame(cbind(va, z, y))
        names(data) = c("v1", "v2", "z", "y")
        
        cohort.data = data[1: n1, ]
        
        data = data[-c(1:n1), ] 
        
        # case control study
        control = data %>%
                filter(y == 0)
        control = control[c(1: (n2 - n21)), ]
        
        case = data %>%
                filter(y == 1)
        
        while (nrow(case) < n21) {
                
                v1 = rep(1, 10000)
                v2 = runif(10000, -2, 2)
                va = cbind(v1, v2)
                
                pz = exp(va %*% eta.true) / (1 + exp(va %*% eta.true))
                z = rbinom(10000, 1, pz)
                
                P.mat = t(mapply(getProbScalarRR, va %*% gamma.true, va %*% beta.true))
                
                y = rep(0, 10000)
                
                for(i in 1:2){
                        y[z == i-1] = rbinom(length(which(z==i-1)), 1, P.mat[z==i-1,i])
                }
                
                data = as.data.frame(cbind(va, z, y))
                names(data) = c("v1", "v2", "z", "y")
                
                case = data %>%
                        filter(y == 1) %>%
                        rbind(case)
        }
        case = case[1: n21, ]
        
        casecrt.dat = rbind(case, control)
        
        return(list(cohort = cohort.data, case.control = casecrt.dat))
     
}

do.one.simulation = function(n1, n2, n21, gamma.true, beta.true, eta.true) {
     #------------- meta dataset generation ------------
     metadata = simData(n1, n2, n21, gamma.true, beta.true, eta.true)
     cohort.dat = metadata$cohort
     casecrt.dat = metadata$case.control
     print(sum(cohort.dat$y)/n1)
     #save(cohort.dat, casecrt.dat, file = paste("data/crt", n1, "_casecrt", n2, "_seed", seed, "_",st, ".RData", sep = ""))
     
     # -------- estimation, corhort only ----------
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
     
     test2 = max.likelihood.casecrt(y = casecrt.dat$y, 
                                    z = casecrt.dat$z, 
                                    va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                    vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                    ps.mat = ps.mat,
                                    alpha.start = c(0,0),
                                    beta.start = c(0,0))
     
     test3 = max.likelihood.casecrt(y = casecrt.dat$y, 
                                    z = casecrt.dat$z, 
                                    va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                    vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                    ps.mat = ps.mat,
                                    alpha.start = test1$point.est[c(1,2)],
                                    beta.start = test1$point.est[c(3,4)])
     # -------- estimation, both ----------
     ps.model1 = glm(c(cohort.dat$z, casecrt.dat$z) ~ c(cohort.dat$v2, casecrt.dat$v2), family = binomial(link='logit'))
     
     ps = exp(predict(ps.model1)) / (exp(predict(ps.model1)) + 1)
     ps = ps[-c(1:nrow(cohort.dat))]
     ps.mat1 = cbind(1-ps, ps)
     
     test4 = max.likelihood.casecrt(y = casecrt.dat$y, 
                                    z = casecrt.dat$z, 
                                    va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                    vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                    ps.mat = ps.mat1,
                                    alpha.start = test1$point.est[c(1,2)],
                                    beta.start = test1$point.est[c(3,4)])
     
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
                                   ps.mat = ps.mat1,
                                   alpha.start = c(0,0),
                                   beta.start = c(0,0))
     
     opt = c(c(test1$coefficients[,c(1,2)]),
             c(test2$coefficients[,c(1,2)]),
             c(test3$coefficients[,c(1,2)]),
             c(test4$coefficients[,c(1,2)]),
             c(test.all$coefficients[,c(1,2)]))
     
     return(opt)
}

N_sim = 100
est.mat = replicate(N_sim, do.one.simulation(n1, n2, n21, gamma.true, beta.true, eta.true))

save(est.mat, file = paste("data/crt", n1, "_casecrt", n2, "_seed", seed, "_", st, ".RData", sep = ""))



