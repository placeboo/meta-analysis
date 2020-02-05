rm(list = ls())  
source("Rcodes/0_library.R")

seed = 2
set.seed(seed)
# ------ study setting ------------
## cohort study
m1 = 1
n1 = 1000

## case-control study
m2 = 1
n2 = 500
n21 = round(n2/3)

# ------- model setting ---------
beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
#beta.true = c(-0.5, 1)
eta.true = c(0.2,-1)
# beta.true = c(-1, 0.5)
# gamma.true = c(0,1)
# eta.true = c(1,-1)

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

#------------- meta dataset generation ------------
metadata = simData(n1, n2, n21, gamma.true, beta.true, eta.true)
cohort.dat = metadata$cohort
casecrt.dat = metadata$case.control
sum(cohort.dat$y)/n1
save(cohort.dat, casecrt.dat, file = paste("data/crt", n1, "_casecrt", n2, "_seed", seed, ".RData", sep = ""))


# test if the simultion is correct for case control study
sum(casecrt.dat$y) / nrow(casecrt.dat) # 0.4


