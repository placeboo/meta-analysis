rm(list = ls())  
source("Rcodes/0_library.R")

seed = 2
set.seed(seed)
# ------ study setting ------------
n = 10000
## cohort study
m1 = 1
n1 = 2000

## case-control study
m2 = 1
n2 = 500
n21 = round(n2/3)

# ------- model setting ---------
beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
#beta.true = c(-0.5, 1)
eta.true = c(0.2,-1)

# ------- data simulation function------
simData = function(n, n1, n2, n21, gamma.true, beta.true, eta.true) {
        
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
        case = case[sample(1:nrow(case), size = n21, replace = TRUE), ]
        
        control = data %>%
                filter(y==0)
        control = control[sample(1:nrow(control), size = n2 - n21, replace = TRUE), ]
        
        casecrt.dat = rbind(case, control)
        
        return(list(cohort = cohort.data, case.control = casecrt.dat))
        
}

#------------- meta dataset generation ------------
metadata = simData(n, n1, n2, n21, gamma.true, beta.true, eta.true)
cohort.dat = metadata$cohort
casecrt.dat = metadata$case.control
sum(cohort.dat$y)/n1
save(cohort.dat, casecrt.dat, file = paste("data/crt", n1, "_casecrt", n2, "_seed", seed, ".RData", sep = ""))


# test if the simultion is correct for case control study
sum(casecrt.dat$y) / nrow(casecrt.dat) # 0.4

