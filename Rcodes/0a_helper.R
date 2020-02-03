# ----- monotonic treatment estimator ------------
max.likelihood.casecrt = function(y, z, va, vb,
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
        cov = matrix(NA, pa+pb+pc, pa+pb+pc)
    } else cov = solve(hessian.mat)
    
    point.est = c(alpha, beta, eta)
    
    opt = orgEst(point.est = point.est,
                 cov = cov,
                 type = "monotone",
                 name = c(paste("alpha", 1:ncol(va), sep = ""), paste("beta", 1:ncol(vb), sep = ""), paste("eta", 1:ncol(va), sep = "")),
                 va=va,
                 vb=vb,
                 coverged = step < max.step)
    
    return(opt)
}

max.likelihood.all= function(case.control = list(y, z, va, vb),
                             cohort = list(y, z, va, vb),
                             alpha.start, beta.start, eta.start,
                             max.step=10^4, thres=10^-6){
    
    pa = length(alpha.start)
    pb = length(beta.start)
    pc = length(eta.start)
    
    casecrt.ny = length(case.control$y)
    cohort.ny = length(cohort$y)
    
    z.uniq = sort(unique(case.control$z))
    nz = length(z.uniq)
    
    
    neg.log.value = function(alpha, beta, eta){
        
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
        casecrt.ps =  exp(case.control$va %*% eta) / (1 + exp(case.control$va %*% eta))
        ps.mat = cbind(1-casecrt.ps, casecrt.ps)
        
        casecrt.P.mat = t(mapply(getProbScalarRR, case.control$va %*% alpha, case.control$vb %*% beta))
        
        casecrt.n_0.vec = apply(casecrt.P.mat, 2, function(x) sum(x==0))
        
        if(identical(all.equal(casecrt.n_0.vec, rep(0, nz)), TRUE)){
            
            py.mat = case.control$y * casecrt.P.mat + (1 - case.control$y)  * (1 - casecrt.P.mat)
            
            for(i in 1: length(z.uniq)){
                value = value - sum(log(py.mat[case.control$z == z.uniq[i],i] * ps.mat[case.control$z == z.uniq[i],i]) - log(apply(py.mat[case.control$z == z.uniq[i], ] * ps.mat[case.control$z == z.uniq[i], ], 1, sum)))
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
        cov = matrix(NA, pa+pb+pc, pa+pb+pc)
    } else cov = solve(hessian.mat)
    
    point.est = c(alpha, beta, eta)
    
    opt = orgEst(point.est = point.est,
                 cov = cov,
                 type = "monotone",
                 name = c(paste("alpha", 1:pa, sep = ""), paste("beta", 1:pb, sep = ""), paste("eta", 1:pc, "")),
                 va=rbind(case.control$va, cohort$va),
                 vb=rbind(case.control$vb,cohort$va),
                 coverged = step < max.step)
    return(opt)
}


orgEst = function(point.est, cov, type, name, va, vb, coverged){
        
        if (any(is.na(cov)) | any(cov < 0)) {
                se.est = conf.lower = conf.upper = p.value = rep(NA, length(point.est))
                
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

