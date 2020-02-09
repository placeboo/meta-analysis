simCohort = function(n, gamma.true, beta.true, eta.true) {
     
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

