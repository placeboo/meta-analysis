#!/usr/local/bin/Rscript
rm(list = ls()) 


library(brm)
library(dplyr)

source("Rcodes/0b_helper_1Dim.R") 


# ------ study setting ------------
gamma.true = c(0.5)
beta.true = c(-7)
eta.true = c(0.2)
n1 = 500
n2 = 200

# ------- data simulation function------
set.seed(1)
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

cohort.dat = simCohort(n1, gamma.true, beta.true, eta.true)

test1 = mem.mono(y = cohort.dat$y,
                 z = cohort.dat$z,
                 va = as.matrix(cohort.dat$va),
                 vb = as.matrix(cohort.dat$va),
                 alpha.start = 0,
                 beta.start = 0)


casectrl.dat = simCasecrt(n2, pr=1/3, gamma.true, beta.true, eta.true)   
ps = mean(cohort.dat$z)
eta = log(ps/(1-ps))


test.all = max.likelihood.all(case.control = list(y = casectrl.dat$y, 
                                                  z = casectrl.dat$z, 
                                                  va = as.matrix(casectrl.dat$va),
                                                  vb = as.matrix(casectrl.dat$va)),
                              cohort = list(y = cohort.dat$y, 
                                            z = cohort.dat$z, 
                                            va = as.matrix(cohort.dat$va),
                                            vb = as.matrix(cohort.dat$va)),
                              alpha.start = test1$point.est[1],
                              beta.start = test1$point.est[2],
                              eta.start = ps)

test.casectrl = max.likelihood.casecrt(y = casectrl.dat$y, 
                                       z = casectrl.dat$z, 
                                       va = as.matrix(casectrl.dat$va),
                                       vb = as.matrix(casectrl.dat$va),
                                       eta = eta,
                                       alpha.start = c(0), 
                                       beta.start = c(0))
# var--------------
var.1d = function(gamma, beta, eta, casectrl.dat, cohort.dat) {
     cohort.va = as.matrix(cohort.dat$va)
     cohort.vb = cohort.va
     cohort.vc = cohort.va
     cohort.Pmat = t(mapply(getProbScalarRR, cohort.va %*% gamma, cohort.vb %*% beta))
     cohort.pi = exp(cohort.vc %*% eta) / (exp(cohort.vc %*% eta) + 1)    
     cohort.p.vec = rep(NA, nrow(cohort.va))
     for (j in 1:nrow(cohort.dat)){
          cohort.p.vec[j] = cohort.Pmat[j, cohort.dat$z[j]+1]
     }
     
     pl1.pgamma = (cohort.dat$y - cohort.p.vec) / (1 - cohort.p.vec) * 
          (cohort.dat$z - (1-cohort.Pmat[,1])/(1-cohort.Pmat[,1] + 1-cohort.Pmat[,2])) * cohort.va
     
     pl1.pbeta =  (cohort.dat$y - cohort.p.vec) / (1 - cohort.p.vec) * ((1-cohort.Pmat[,1]) * (1 - cohort.Pmat[,2])) / ((1-cohort.Pmat[,1]) + (1 - cohort.Pmat[,2])) * cohort.vb
     
     pl1.peta = 0 
     
     # l2 deviatives -------------
     pl2.pgamma = 0
     pl2.pbeta = 0
     pl2.peta = (cohort.dat$z - cohort.pi) * cohort.vc 
     
     # l3 deviatives ------------
     casectrl.va = as.matrix(casectrl.dat$va)
     casectrl.vb = casectrl.va
     casectrl.vc = casectrl.va
     casectrl.Pmat = t(mapply(getProbScalarRR, casectrl.va %*% gamma, casectrl.vb %*% beta))
     casectrl.pi = exp(casectrl.vc %*% eta) / (exp(casectrl.vc %*% eta) + 1)
     
     casectrl.p.vec = rep(NA, nrow(casectrl.dat))
     for (j in 1:nrow(casectrl.dat)){
          casectrl.p.vec[j] = casectrl.Pmat[j, casectrl.dat$z[j]+1]
     }
     qA = (casectrl.dat$y *  casectrl.Pmat[,2] + (1 - casectrl.dat$y) * (1 -  casectrl.Pmat[,2])) * casectrl.pi
     qB = qA + (casectrl.dat$y *  casectrl.Pmat[,1] + (1 - casectrl.dat$y) * (1 -  casectrl.Pmat[,1])) * (1 - casectrl.pi)
     casectrl.q = qA/qB 
     
     tmp = (casectrl.dat$z - casectrl.q) / (casectrl.q * (1 -casectrl.q))
     
     pp1.pgamma = casectrl.Pmat[,2] * (1 - casectrl.Pmat[,2])/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2])) * casectrl.va
     pp0.pgamma = casectrl.Pmat[,1] * ( - (1 - casectrl.Pmat[,1])/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.va
     
     pp1.pbeta = casectrl.Pmat[,2] * (((1 - casectrl.Pmat[,1]) * (1 - casectrl.Pmat[,2]))/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.vb
     pp0.pbeta = casectrl.Pmat[,1] * (((1 - casectrl.Pmat[,1]) * (1 - casectrl.Pmat[,2]))/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.vb
     
     pA.pgamma = (2 * casectrl.dat$y - 1)  * pp1.pgamma * casectrl.pi
     pA.pbeta = (2 * casectrl.dat$y - 1)  * pp1.pbeta * casectrl.pi
     pA.peta = (casectrl.dat$y * casectrl.Pmat[,2] + (1 - casectrl.dat$y) * (1 - casectrl.Pmat[,2])) * casectrl.pi * (1 - casectrl.pi) * casectrl.vc
     
     pB.pgamma = pA.pgamma + (2 * casectrl.dat$y - 1) * pp0.pgamma * (1-casectrl.pi)
     pB.pbeta = pA.pbeta + (2 * casectrl.dat$y - 1) * pp0.pbeta * (1-casectrl.pi)
     pB.peta = pA.peta -  (casectrl.dat$y * casectrl.Pmat[,1] + (1 - casectrl.dat$y) * (1 - casectrl.Pmat[,1])) *  casectrl.pi * (1 - casectrl.pi) * casectrl.vc
     
     pq.pgamma = (qB * pA.pgamma - qA * pB.pgamma) * 1/qB * 1/qB
     pq.pbeta = (qB * pA.pbeta - qA * pB.pbeta) * 1/qB * 1/qB
     pq.peta = (qB * pA.peta- qA * pB.peta) * 1/qB * 1/qB
     
     pl3.pgamma = tmp * pq.pgamma
     pl3.pbeta = tmp * pq.pbeta
     pl3.peta = tmp * pq.peta
     
     # combina
     dl.dgamma = c(pl1.pgamma, pl3.pgamma)
     dl.dbeta = c(pl1.pbeta, pl3.pbeta)
     dl.deta = c(pl2.peta, pl3.peta)
     
     dl = cbind(dl.dgamma, dl.dbeta, dl.deta)
     hessian.mat = t(dl) %*% dl
     var.mat = solve(hessian.mat)
     return(var.mat)
}

var.casectrl.1d = function(gamma, beta, eta, casectrl.dat) {
     # l3 deviatives ------------
     casectrl.va = as.matrix(casectrl.dat$va)
     casectrl.vb = casectrl.va
     casectrl.vc = casectrl.va
     casectrl.Pmat = t(mapply(getProbScalarRR, casectrl.va %*% gamma, casectrl.vb %*% beta))
     casectrl.pi = exp(casectrl.vc %*% eta) / (exp(casectrl.vc %*% eta) + 1)
     
     casectrl.p.vec = rep(NA, nrow(casectrl.dat))
     for (j in 1:nrow(casectrl.dat)){
          casectrl.p.vec[j] = casectrl.Pmat[j, casectrl.dat$z[j]+1]
     }
     qA = (casectrl.dat$y *  casectrl.Pmat[,2] + (1 - casectrl.dat$y) * (1 -  casectrl.Pmat[,2])) * casectrl.pi
     qB = qA + (casectrl.dat$y *  casectrl.Pmat[,1] + (1 - casectrl.dat$y) * (1 -  casectrl.Pmat[,1])) * (1 - casectrl.pi)
     casectrl.q = qA/qB 
     
     tmp = (casectrl.dat$z - casectrl.q) / (casectrl.q * (1 -casectrl.q))
     
     pp1.pgamma = casectrl.Pmat[,2] * (1 - casectrl.Pmat[,2])/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2])) * casectrl.va
     pp0.pgamma = casectrl.Pmat[,1] * ( - (1 - casectrl.Pmat[,1])/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.va
     
     pp1.pbeta = casectrl.Pmat[,2] * (((1 - casectrl.Pmat[,1]) * (1 - casectrl.Pmat[,2]))/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.vb
     pp0.pbeta = casectrl.Pmat[,1] * (((1 - casectrl.Pmat[,1]) * (1 - casectrl.Pmat[,2]))/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.vb
     
     pA.pgamma = (2 * casectrl.dat$y - 1)  * pp1.pgamma * casectrl.pi
     pA.pbeta = (2 * casectrl.dat$y - 1)  * pp1.pbeta * casectrl.pi
     pA.peta = (casectrl.dat$y * casectrl.Pmat[,2] + (1 - casectrl.dat$y) * (1 - casectrl.Pmat[,2])) * casectrl.pi * (1 - casectrl.pi) * casectrl.vc
     
     pB.pgamma = pA.pgamma + (2 * casectrl.dat$y - 1) * pp0.pgamma * (1-casectrl.pi)
     pB.pbeta = pA.pbeta + (2 * casectrl.dat$y - 1) * pp0.pbeta * (1-casectrl.pi)
     pB.peta = pA.peta -  (casectrl.dat$y * casectrl.Pmat[,1] + (1 - casectrl.dat$y) * (1 - casectrl.Pmat[,1])) *  casectrl.pi * (1 - casectrl.pi) * casectrl.vc
     
     pq.pgamma = (qB * pA.pgamma - qA * pB.pgamma) * 1/qB * 1/qB
     pq.pbeta = (qB * pA.pbeta - qA * pB.pbeta) * 1/qB * 1/qB
     pq.peta = (qB * pA.peta- qA * pB.peta) * 1/qB * 1/qB
     
     pl3.pgamma = tmp * pq.pgamma
     pl3.pbeta = tmp * pq.pbeta
     pl3.peta = tmp * pq.peta
     
     # combina

     dl = cbind(pl3.pgamma, pl3.pbeta, pl3.peta)
     hessian.mat = t(dl) %*% dl
     var.mat = solve(hessian.mat)
     return(var.mat)
}

var.mat = var.1d(gamma = test.all$point.est[1],
       beta = test.all$point.est[2],
       eta = test.all$point.est[3],
       casectrl.dat,
       cohort.dat)

test.all$cov

sqrt(diag(var.mat))/test.all$se.est

var.casectrl.1d(gamma = test.casectrl$point.est[1],
                beta = test.casectrl$point.est[2],
                eta = eta,
                casectrl.dat)

gamma = test.all$point.est[1]
beta = test.all$point.est[2]
eta = test.all$point.est[3]

# cohort deviatives
cohort.va = as.matrix(cohort.dat$va)
cohort.vb = cohort.va
cohort.vc = cohort.va
cohort.Pmat = t(mapply(getProbScalarRR, cohort.va %*% gamma, cohort.vb %*% beta))
cohort.pi = exp(cohort.vc %*% eta) / (exp(cohort.vc %*% eta) + 1)

cohort.p.vec = rep(NA, nrow(cohort.dat))
for (j in 1:nrow(cohort.dat)){
     cohort.p.vec[j] = cohort.Pmat[j, cohort.dat$z[j]+1]
}

pl1.pgamma = (cohort.dat$y - cohort.p.vec) / (1 - cohort.p.vec) * 
                        (cohort.dat$z - (1-cohort.Pmat[,1])/(1-cohort.Pmat[,1] + 1-cohort.Pmat[,2])) * cohort.va

pl1.pbeta =  (cohort.dat$y - cohort.p.vec) / (1 - cohort.p.vec) * ((1-cohort.Pmat[,1]) * (1 - cohort.Pmat[,2])) / ((1-cohort.Pmat[,1]) + (1 - cohort.Pmat[,2])) * cohort.vb

pl1.peta = 0 

# l2 deviatives -------------
pl2.pgamma = 0
pl2.pbeta = 0
pl2.peta = (cohort.dat$z - cohort.pi) * cohort.vc 

# l3 deviatives ------------
casectrl.va = as.matrix(casectrl.dat$va)
casectrl.vb = casectrl.va
casectrl.vc = casectrl.va
casectrl.Pmat = t(mapply(getProbScalarRR, casectrl.va %*% gamma, casectrl.vb %*% beta))
casectrl.pi = exp(casectrl.vc %*% eta) / (exp(casectrl.vc %*% eta) + 1)

casectrl.p.vec = rep(NA, nrow(casectrl.dat))
for (j in 1:nrow(casectrl.dat)){
     casectrl.p.vec[j] = casectrl.Pmat[j, casectrl.dat$z[j]+1]
}
qA = (casectrl.dat$y *  casectrl.Pmat[,2] + (1 - casectrl.dat$y) * (1 -  casectrl.Pmat[,2])) * casectrl.pi
qB = qA + (casectrl.dat$y *  casectrl.Pmat[,1] + (1 - casectrl.dat$y) * (1 -  casectrl.Pmat[,1])) * (1 - casectrl.pi)
casectrl.q = qA/qB 

tmp = (casectrl.dat$z - casectrl.q) / (casectrl.q * (1 -casectrl.q))

pp1.pgamma = casectrl.Pmat[,2] * (1 - casectrl.Pmat[,2])/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2])) * casectrl.va
pp0.pgamma = casectrl.Pmat[,1] * ( - (1 - casectrl.Pmat[,1])/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.va

pp1.pbeta = casectrl.Pmat[,2] * (((1 - casectrl.Pmat[,1]) * (1 - casectrl.Pmat[,2]))/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.vb
pp0.pbeta = casectrl.Pmat[,1] * (((1 - casectrl.Pmat[,1]) * (1 - casectrl.Pmat[,2]))/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.vb

pA.pgamma = (2 * casectrl.dat$y - 1)  * pp1.pgamma * casectrl.pi
pA.pbeta = (2 * casectrl.dat$y - 1)  * pp1.pbeta * casectrl.pi
pA.peta = (casectrl.dat$y * casectrl.Pmat[,2] + (1 - casectrl.dat$y) * (1 - casectrl.Pmat[,2])) * casectrl.pi * (1 - casectrl.pi) * casectrl.vc

pB.pgamma = pA.pgamma + (2 * casectrl.dat$y - 1) * pp0.pgamma * (1-casectrl.pi)
pB.pbeta = pA.pbeta + (2 * casectrl.dat$y - 1) * pp0.pbeta * (1-casectrl.pi)
pB.peta = pA.peta -  (casectrl.dat$y * casectrl.Pmat[,1] + (1 - casectrl.dat$y) * (1 - casectrl.Pmat[,1])) *  casectrl.pi * (1 - casectrl.pi) * casectrl.vc

pq.pgamma = (qB * pA.pgamma - qA * pB.pgamma) * 1/qB * 1/qB
pq.pbeta = (qB * pA.pbeta - qA * pB.pbeta) * 1/qB * 1/qB
pq.peta = (qB * pA.peta- qA * pB.peta) * 1/qB * 1/qB

pl3.pgamma = tmp * pq.pgamma
pl3.pbeta = tmp * pq.pbeta
pl3.peta = tmp * pq.peta

# combina
dl.dgamma = c(pl1.pgamma, pl3.pgamma)
dl.dbeta = c(pl1.pbeta, pl3.pbeta)
dl.deta = c(pl2.peta, pl3.peta)

dl = cbind(dl.dgamma, dl.dbeta, dl.deta)
hessian.mat = t(dl) %*% dl
var.mat = solve(hessian.mat)
diag(hessian.mat)/diag(solve(test.all$cov))
sqrt(diag(var.mat))/test.all$se.est

test.all$cov
var.mat

test.casectrl$cov



# check with 
gamma.test = test1$point.est[1]
beta.test = test1$point.est[2]

cohort.Pmat = t(mapply(getProbScalarRR, cohort.va %*% gamma.test, cohort.vb %*% beta.test))

cohort.p.vec = rep(NA, nrow(cohort.dat))
for (j in 1:nrow(cohort.dat)){
     cohort.p.vec[j] = cohort.Pmat[j, cohort.dat$z[j]+1]
}

pl1.pgamma = (cohort.dat$y - cohort.p.vec) / (1 - cohort.p.vec) * 
     (cohort.dat$z - (1-cohort.Pmat[,1])/(1-cohort.Pmat[,1] + 1-cohort.Pmat[,2])) * cohort.va

pl1.pbeta =  (cohort.dat$y - cohort.p.vec) / (1 - cohort.p.vec) * ((1-cohort.Pmat[,1]) * (1 - cohort.Pmat[,2])) / ((1-cohort.Pmat[,1]) + (1 - cohort.Pmat[,2])) * cohort.vb

t(cbind(pl1.pgamma, pl1.pbeta)) %*% cbind(pl1.pgamma, pl1.pbeta)
solve(test1$cov)

