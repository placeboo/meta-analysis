rm(list = ls())  
source("Rcodes/0_library.R")

seed = 1
set(seed)
# ------ study setting ------------
n = 10000
## cohort study
m1 = 1
n1 = 1000

## case-control study
m2 = 1
n2 = 500
n20 = 200

# ------- model setting ---------
gamma.true = c(0, 0.2)
beta.true = c(-5, 0.5)
eta.true = c(0.2,-1)

# ------- data simulation function------
simData = function(n, n1, n2, n20, gamma.true, beta.true, eta.true) {
        
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
        
        
        # cohort data 
        cohort.idx = sample(1: n, size = n1, replace = TRUE)
        cohort.data = data[cohort.idx,]
        
        # case control study
        case = data %>%
                filter(y==1)
        case = case[sample(1:nrow(case), size = n2-n20, replace = TRUE), ]
        
        control = data %>%
                filter(y==0)
        control = control[sample(1:nrow(control), size = n20, replace = TRUE), ]
        
        casecrt.dat = rbind(case, control)
        
        return(list(cohort = cohort.data, case.control = casecrt.dat))
        
}

#------------- meta dataset generation ------------
metadata = simData(n, n1, n2, n20, gamma.true, beta.true, eta.true)
cohort.dat = metadata$cohort
casecrt.dat = metadata$case.control

save(cohort.dat, casecrt.dat, file = paste("data/meta-data-generation-seed", seed, ".RData", sep = ""))


