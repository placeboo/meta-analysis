#!/usr/local/bin/Rscript
rm(list = ls())

args = commandArgs(TRUE)

n = as.numeric(args[[1]])
seed = as.numeric(args[[2]])
st = "01302020"
library(dplyr)
library(brm)

orgEst = function(point.est, cov, type, name, va, vb, coverged) {
     
     if (any(is.na(cov))) {
          se.est = conf.lower = conf.upper = p.value = rep(NA, ncol(va) + ncol(vb))
          
     } else {
          se.est = sqrt(diag(cov))
          conf.lower = point.est + stats::qnorm(0.025) * se.est
          conf.upper = point.est + stats::qnorm(0.975) * se.est
          
          p.temp = stats::pnorm(point.est/se.est, 0, 1)
          p.value = 2 * pmin(p.temp, 1 - p.temp)
     }
     
     
     if (type == "monotone") {
          names(point.est) = names(se.est) = rownames(cov) = colnames(cov) = names(conf.lower) = names(conf.upper) = names(p.value) = name
          coefficients = cbind(point.est, se.est, conf.lower, conf.upper,
                               p.value)
          linear.predictors = va %*% point.est[1:ncol(va)]
          param.est = exp(linear.predictors)
          
          output = list(point.est = point.est,
                        se.est = se.est,
                        cov = cov,
                        conf.lower = conf.lower,
                        conf.upper = conf.upper,
                        p.value = p.value,
                        coefficients = coefficients,
                        param = param.est,
                        converged = coverged)
          
          class(output) = c("mem", type, "list")
          
     }
     return(output)
}


simData = function(n, n1, gamma.true, beta.true, eta.true) {
     
     v1 = rep(1, 20000)
     v2 = runif(20000, -2, 2)
     va = cbind(v1, v2)
     
    
     z = rbinom(20000, 1, 0.5)
     
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
          
          z = rbinom(10000, 1, 0.5)
          
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
max.likelihood.casecrt = function(y, z, va, vb,
                                  alpha.start, beta.start,
                                  max.step=10^4, thres=10^-6) {
     if (!is.numeric(z)) stop("'z' should be numerical!")
     
     pa = length(alpha.start)
     pb = length(beta.start)
     ny = length(y)
     
     z.uniq = sort(unique(z))
     nz = length(z.uniq)
     
     neg.log.value = function(alpha, beta){
          
          Pzmin_Pzmax = t(mapply(getProbScalarRR, va %*% alpha, vb %*% beta))
          P.mat = matrix(0, ncol = nz, nrow = ny)
          P.mat[, c(1, nz)] = Pzmin_Pzmax
          
          P.mat[, -c(1, nz)] = Pzmin_Pzmax[,1] * exp(va %*% alpha %*% t(z.uniq)[-c(1,nz)])
          
          n_0.vec = apply(P.mat, 2, function(x) sum(x==0))
          
          if(identical(all.equal(n_0.vec, rep(0, nz)), TRUE)){
               value = 0
               py.mat = y * P.mat + (1 - y)  * (1 - P.mat)
               
               for (i in 1: length(z.uniq)){
                    value = value - sum(log(py.mat[z == z.uniq[i],i]) - log(apply(py.mat[z == z.uniq[i], ], 1, sum)))
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
          
          return(neg.log.value(alpha, beta))
     }
     
     neg.log.likelihood.alpha = function(alpha){
          return(neg.log.value(alpha, beta))
     }
     
     neg.log.likelihood.beta = function(beta){
          
          return(neg.log.value(alpha, beta))
     }
     
     Diff = function(x, y) sum((x - y)^2)/sum(x^2 + thres)
     alpha = alpha.start
     beta = beta.start
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
          
          diff = max(diff1, diff2)
     }
     # have hessian matrix
     hessian.mat = try(optimHess(c(alpha,beta), neg.log.likelihood), silent = TRUE)
     
     if ("try-error" %in% class(hessian.mat)) {
          return( alpha.cov = matrix(NA, pa+pb, pa+pb))
     }
     
     cov = solve(hessian.mat)
     point.est = c(alpha, beta)
     
     opt = orgEst(point.est = point.est,
                  cov = cov,
                  type = "monotone",
                  name = c(paste("alpha", 1:ncol(va), sep = ""), paste("beta", 1:ncol(vb), sep = "")),
                  va=va,
                  vb=vb,
                  coverged = step < max.step)
     
     return(opt)
}




# ----- data simulation ----------------
beta.true = c(-7, 0.5) 
gamma.true = c(0,-1)
n1 = round(n/3)

do.one.simulation = function(n, n1, gamma.true, beta.true){
     casecrt.dat = simData(n, n1, gamma.true, beta.true)
     test = max.likelihood.casecrt(y = casecrt.dat$y, 
                                   z = casecrt.dat$z, 
                                   va = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                   vb = cbind(casecrt.dat$v1, casecrt.dat$v2),
                                   
                                   alpha.start = c(0,0),
                                   beta.start = c(-7,0))
     return(c(test$coefficients[,c(1,2)]))
}

N_sim = 100
est.mat = replicate(N_sim, do.one.simulation(n, n1, gamma.true, beta.true))

save(est.mat, file = paste("case-control/casecontrol", n, "_seed", seed, "_", st, ".RData", sep = ""))


