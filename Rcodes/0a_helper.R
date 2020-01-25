# ----- monotonic treatment estimator ------------
max.likelihood.casecrt = function(y, z, va, vb, ps.mat,
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
                        
                        for(i in 1: length(z.uniq)){
                                for (j in c(0,1)) {
                                        idx = (z == z.uniq[i]) & (y == j)
                                        ps = ps.mat[idx, ]
                                        py = j * P.mat[idx, ] + (1 - j) * (1 - P.mat[idx, ])
                                        value = value - sum(log(py[, i]) + log(ps[,i]) - log(apply(ps * py, 1, sum)))
                                }
       
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
                return()
                alpha.cov = matrix(NA, pa+pb, pa+pb)
        }
        
        cov = solve(hessian.mat)
        point.est = c(alpha, beta)
        
        likelihood = exp(-neg.log.likelihood(point.est))
        opt = orgEst(point.est = point.est,
                     cov = cov,
                     type = "monotone",
                     name = c(paste("alpha", 1:ncol(va), sep = ""), paste("beta", 1:ncol(vb), sep = "")),
                     va=va,
                     vb=vb,
                     coverged = step < max.step,
                     likelihood = likelihood)
        
        return(opt)
}

max.likelihood.all= function(case.control = list(y, z, va, vb),
                             cohort = list(y, z, va, vb),
                             ps.mat,
                             alpha.start, beta.start, max.step=10^5, thres=10^-5){
     
     
    
     pa = length(alpha.start)
     pb = length(beta.start)
     
     casecrt.ny = length(case.control$y)
     cohort.ny = length(cohort$y)
     
     z.uniq = sort(unique(case.control$z))
     nz = length(z.uniq)
     
     
     neg.log.value = function(alpha, beta){
             
             # cohort
             cohort.Pzmin_Pzmax = t(mapply(getProbScalarRR, cohort$va %*% alpha, cohort$vb %*% beta))
             cohort.P.mat = matrix(0, ncol = nz, nrow = cohort.ny)
             cohort.P.mat[, c(1, nz)] = cohort.Pzmin_Pzmax
             
             cohort.P.mat[, -c(1, nz)] = cohort.Pzmin_Pzmax[,1] * exp(cohort$va %*% alpha %*% t(z.uniq)[-c(1,nz)])
             
             cohort.n_0.vec = apply(cohort.P.mat, 2, function(x) sum(x==0))
             
             if(identical(all.equal(cohort.n_0.vec, rep(0, nz)), TRUE)){
                     value = 0
                     
                     for(i in 1: length(z.uniq)){
                             value = value - sum(cohort$y[cohort$z==z.uniq[i]] * log(cohort.P.mat[cohort$z==z.uniq[i],i]) + (1-cohort$y[cohort$z==z.uniq[i]]) * log(1-cohort.P.mat[cohort$z==z.uniq[i],i]))
                             
                     }
             }
             else{
                     return(10000)
             }
             
             # case control
             casecrt.Pzmin_Pzmax = t(mapply(getProbScalarRR, case.control$va %*% alpha, case.control$vb %*% beta))
             
             casecrt.P.mat = matrix(0, ncol = nz, nrow = casecrt.ny)
             casecrt.P.mat[, c(1, nz)] = casecrt.Pzmin_Pzmax
             
             casecrt.P.mat[, -c(1, nz)] = casecrt.Pzmin_Pzmax[,1] * exp(case.control$va %*% alpha %*% t(z.uniq)[-c(1,nz)])
             
             casecrt.n_0.vec = apply(casecrt.P.mat, 2, function(x) sum(x==0))
             
             if(identical(all.equal(casecrt.n_0.vec, rep(0, nz)), TRUE)){
                
                     for(i in 1: length(z.uniq)){
                             for (j in c(0,1)) {
                                     idx = (case.control$z == z.uniq[i]) & (case.control$y == j)
                                     ps = ps.mat[idx, ]
                                     py = j * casecrt.P.mat[idx, ] + (1 - j) * (1 - casecrt.P.mat[idx, ])
                                     value = value - sum(log(py[, i]) + log(ps[,i]) - log(apply(ps * py, 1, sum)))
                             }
                             
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
             return()
             alpha.cov = matrix(NA, pa+pb, pa+pb)
     }
     
     cov = solve(hessian.mat)
     point.est = c(alpha, beta)
     
     likelihood = exp(-neg.log.likelihood(point.est))
     opt = orgEst(point.est = point.est,
                  cov = cov,
                  type = "monotone",
                  name = c(paste("alpha", 1:ncol(case.control$va), sep = ""), paste("beta", 1:ncol(case.control$vb), sep = "")),
                  va=case.control$va,
                  vb=case.control$vb,
                  coverged = step < max.step,
                  likelihood = likelihood)
     
     return(opt)
}

orgEst = function(point.est, cov, type, name, va, vb, coverged, likelihood) {
        
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
                              converged = coverged,
                              likelihood = likelihood)
                
                class(output) = c("mem", type, "list")
                
        }
        return(output)
}

