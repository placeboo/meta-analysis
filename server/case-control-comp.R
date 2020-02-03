rm(list = ls())

library(dplyr)
library(brm)


simData = function(n, n1, gamma.true, beta.true, eta.true) {
     
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
     
     # case control study
     control = data %>%
          filter(y == 0)
     control = control[c(1: (n - n1)), ]
     
     case = data %>%
          filter(y == 1)
     
     while (nrow(case) < n1) {
          
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
     case = case[1: n1, ]
     casecrt.dat = rbind(case, control)
     return(casecrt.dat)
}

# ----- monotonic treatment estimator ------------

max.likelihood.casecrt2 = function(y, z, va, vb,
                                   alpha.start, beta.start, eta.start,
                                   max.step=10^4, thres=10^-6) {
     if (!is.numeric(z)) stop("'z' should be numerical!")
     
     pa = length(alpha.start)
     pb = length(beta.start)
     pc = length(eta.start)
     
     ny = length(y)
     
     z.uniq = sort(unique(z))
     nz = length(z.uniq)
     
     neg.log.value = function(alpha, beta, eta){
          
          ps = exp(va %*% eta) / (1 + exp(va %*% eta))
          ps.mat = cbind(1-ps, ps)
          
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
     
     neg.log.likelihood.eta = function(eta){
          return(neg.log.value(alpha, beta, eta))
     }
     
     neg.log.likelihood.beta = function(beta){
          return(neg.log.value(alpha, beta, eta))
     }
     
     Diff = function(x, y) sum((x - y)^2)/sum(x^2 + thres)
     alpha = alpha.start
     beta = beta.start
     eta = eta.start
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
          
          opt3 = stats::optim(eta, neg.log.likelihood.eta, control = list(maxit = max.step))
          diff3 = Diff(opt3$par, eta)
          eta = opt3$par
          
          diff = max(diff1, diff2, diff3)
     }
     # have hessian matrix
     hessian.mat = try(optimHess(c(alpha,beta,eta), neg.log.likelihood), silent = TRUE)
     
     if ("try-error" %in% class(hessian.mat)) {
          se = c(NA, pa+pb+pc)
     } else se = sqrt(diag(solve(hessian.mat)))
     
     point.est = c(alpha, beta, eta)
     return(c(point.est, se))
}

# ----- data simulation ----------------
do.one.simulation = function(n, n1, gamma.true, beta.true, eta.true){
     casecrt.dat = simData(n, n1, gamma.true, beta.true, eta.true)
     test = max.likelihood.casecrt2(y = casecrt.dat$y, 
                                     z = casecrt.dat$z, 
                                     va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                     vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                     alpha.start = c(0,0),
                                     beta.start = c(0,0),
                                     eta.start = c(0,0))
     
}
n = 5000
seed = 3
set.seed(seed)

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)

n1 = round(n/3)

sim_N = 1
est.mat = replicate(sim_N, do.one.simulation(n, n1, gamma.true, beta.true, eta.true))
