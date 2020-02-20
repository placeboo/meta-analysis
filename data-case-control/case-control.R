#!/usr/local/bin/Rscript
rm(list = ls()) 

args = commandArgs(TRUE)

n = as.numeric(args[[1]])
seed = as.numeric(args[[2]])

library(brm)
library(dplyr)
source("Rcodes/0c_simulation.R")
source("Rcodes/0a_helper.R")

# ------ study setting ------------
st = "02182020"

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)

max.likelihood.casectrl = function(y, z, va, vb,
                                   alpha.start, beta.start, eta.start,
                                   max.step=10^4, thres=10^-6) {
     
     pa = length(alpha.start)
     pb = length(beta.start)
     pc = length(eta.start)
     
     
     ny = length(y)
     
     z.uniq = sort(unique(z))
     nz = length(z.uniq)
     
     eta = eta.start #!!! Trick
     
     neg.log.value = function(alpha, beta, eta){
          
          casecrt.ps =  exp(va %*% eta) / (1 + exp(va %*% eta))
          ps.mat = cbind(1-casecrt.ps, casecrt.ps)
          
          Pzmin_Pzmax = t(mapply(getProbScalarRR, va %*% alpha, vb %*% beta))
          P.mat = matrix(0, ncol = nz, nrow = ny)
          P.mat[, c(1, nz)] = Pzmin_Pzmax
          
          P.mat[, -c(1, nz)] = Pzmin_Pzmax[,1] * exp(va %*% alpha %*% t(z.uniq)[-c(1,nz)])
          
          n_0.vec = apply(P.mat, 2, function(x) sum(x==0))
          
          if(identical(all.equal(n_0.vec, rep(0, nz)), TRUE)){
               value = 0
               py.mat = y * P.mat + (1 - y)  * (1 - P.mat)
               
               for (i in 1: length(z.uniq)){
                    value = value - sum(log(py.mat[z == z.uniq[i],i] * ps.mat[z == z.uniq[i],i]) - log(apply(py.mat[z == z.uniq[i], ] * ps.mat[z == z.uniq[i], ], 1, sum)))
               }
               return(value)
          }
          else{
               return(10000)
          }
     }
     
     
     neg.log.likelihood = function(pars){
          alpha = pars[1:pa]
          beta = pars[(pa + 1):(pa + pb)]
          eta = pars[(pa + pb + 1) : (pa + pb + pc)]
          
          return(neg.log.value(alpha, beta, eta))
     }
     
     neg.log.likelihood.alpha = function(alpha){
          return(neg.log.value(alpha, beta, eta))
     }
     
     # neg.log.likelihood.eta = function(eta){
     #      return(neg.log.value(alpha, beta, eta))
     # }
     
     neg.log.likelihood.beta = function(beta){
          return(neg.log.value(alpha, beta, eta))
     }
     
     Diff = function(x, y) sum((x - y)^2)/sum(x^2 + thres)
     alpha = alpha.start
     beta = beta.start
     # eta = eta.start
     diff = thres + 1
     step = 0
     
     
     while(diff > thres & step < max.step){
          step = step + 1
          opt1 = stats::optim(alpha, neg.log.likelihood.alpha,
                              control = list(maxit = max.step))
          diff1 = Diff(opt1$par, alpha)
          alpha = opt1$par
          
          opt2 = stats::optim(beta, neg.log.likelihood.beta, control = list(maxit = max.step))
          
          diff2 = Diff(opt2$par, beta)
          beta = opt2$par
          
          # opt3 = stats::optim(eta, neg.log.likelihood.eta, control = list(maxit = max.step))
          # diff3 = Diff(opt3$par, eta)
          # eta = opt3$par
          
          # diff = max(diff1, diff2, diff3)
          diff = max(diff1, diff2)
     }
     # have hessian matrix
     hessian.mat = try(optimHess(c(alpha,beta,eta.start), neg.log.likelihood), silent = TRUE)
     
     
     if ("try-error" %in% class(hessian.mat)) {
          cov = matrix(NA, pa+pb+pc, pa+pb+pc)
     } else cov = solve(hessian.mat)
     
     point.est = c(alpha, beta, eta.start)
     
     opt = orgEst(point.est = point.est,
                  cov = cov,
                  type = "monotone",
                  name = c(paste("alpha", 1:pa, sep = ""), paste("beta", 1:pb, sep = ""), paste("eta", 1:pc, "")),
                  va=va,
                  vb=va,
                  coverged = step < max.step)
     
     return(opt)
}

do.one.simulation = function(n, pr, gamma.true, beta.true, eta.true) {
     # cohort have eta
     cohort.dat = simCohort(n = 1000, gamma.true, beta.true, eta.true)
     ps.model = glm(cohort.dat$z ~ cohort.dat$v2, family = binomial(link='logit'))
     
     casecrt.dat = simCasecrt(n, pr, gamma.true, beta.true,  eta.true) 
     
     test = max.likelihood.casectrl(y = casecrt.dat$y, 
                                    z = casecrt.dat$z, 
                                    va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                    vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                    alpha.start = c(0,0),
                                    eta.start =ps.model$coefficients,
                                    beta.start = c(0,0))
     
     opt = c(test$coefficients[, c(1,2)])
     
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


save(est.mat, file = paste("data-case-control/casecrt", n, "_seed", seed, "_", st, ".RData", sep = ""))
 