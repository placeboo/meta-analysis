rm(list = ls())  
source("Rcodes/0_library.R")

seed = 1
set(seed)
# ------ study setting ------------
n.min = 51
n.max = 200

## cohort study
m1 = 15

## case-control study
m2 = 15
p.y = 0.1

# ------- model setting ---------
gamma.true = c(0, 1)
beta.true = c(1, -0.5)
n1 = sample(c(n.min:n.max), m1, replace = TRUE)
n2 = sample(c(n.min:n.max), m2, replace = TRUE)

# ------- data simulation function------

# generate from cohort study
simCht = function(n, gamma.true, beta.true) {

     v1 = rep(1, n)
     v2 = runif(n, -2, 2)
     va = cbind(v1, v2)
     z = sample(x = c(0:2), size = n, replace = T)
     
     
     Pzmin_Pzmax = t(mapply(getProbScalarRR, va %*% gamma.true * 2, va %*% beta.true))
     
     P.mat = matrix(0, ncol = 3, nrow = n)
     P.mat[, c(1, 3)] = Pzmin_Pzmax
     P.mat[, 2] = Pzmin_Pzmax[,1] * exp(va %*% gamma.true)
     
     y = rep(0, n)
     
     for(i in 1:3){
          y[z == i-1] = rbinom(length(which(z==i-1)), 1, P.mat[z==i-1,i])
     }
     
     mat = cbind(va, z, y)
     colnames(mat) = c("v1", "v2", "z", "y")
     return(mat)
}

# case-control study 
simCc = function(p.y, n, gamma.true, beta.true) {
     
     v1 = rep(1, n)
     v2 = runif(n, -2, 2)
     va = cbind(v1, v2)
     
     # generate y
     y = rbinom(n, size = 1, prob = p.y) 
     # At lease one case
     if (sum(y) == 0) {
          y[sample(1:length(y), size = 1)] = 1
     }
     
     Pzmin_Pzmax = t(mapply(getProbScalarRR, va %*% gamma.true * 2, va %*% beta.true))
     
     P.mat = matrix(0, ncol = 3, nrow = n)
     P.mat[, c(1, 3)] = Pzmin_Pzmax
     P.mat[, 2] = Pzmin_Pzmax[,1] * exp(va %*% gamma.true)
     
     # generate Z
     ## when Y = 1
     pz_y1 = apply(P.mat, 1, function(x) x/sum(x)) 
     ## when Y = 0
     pz_y0 = apply(P.mat, 1, function(x) (1 - x)/sum(1-x)) 
     
     z = rep(0, n)
     
     z[y == 1] = apply(as.matrix(pz_y1[ , y == 1]), 2, function(x) sample(x = 0:2, size =1, prob = x))
     z[y == 0] = apply(as.matrix(pz_y0[ , y == 0]), 2, function(x) sample(x = 0:2, size =1, prob = x))
     
     mat = cbind(va, z, y)
     colnames(mat) = c("v1", "v2", "z", "y")
     
     return(mat)
} 

#------------- meta dataset generation ------------
cohort.dat = NULL
casecrt.dat = NULL

for (i in 1: m1) {
     
     tmp = simCht(n1[i], gamma.true, beta.true)
     tmp = as.data.frame(tmp)
     tmp$trial = i
     
     cohort.dat = rbind(cohort.dat, tmp)
}

for (i in 1: m2){
     
     tmp = simCc(p.y, n2[i], gamma.true, beta.true)
     tmp = as.data.frame(tmp)
     tmp$trial = i + m1
     
     casecrt.dat = rbind(casecrt.dat, tmp)
}


cohort.dat %>%
     group_by(trial) %>%
     count()

casecrt.dat %>%
     group_by(trial) %>%
     count()

