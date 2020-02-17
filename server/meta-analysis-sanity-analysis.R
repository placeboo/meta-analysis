rm(list = ls())

st = "02082020"

n1.vec = c(500, 1000, 2000, 4000, 8000)
n2.vec =  c(200, 500, 1000, 2000, 4000)

beta.true = c(-7)
gamma.true = c(0.5)
eta.true = c(0.2)
seed.vec = 1:20

final = matrix(NA, ncol = 3)

for (n1 in n1.vec) {

     mle_tmp_dat = matrix(NA, nrow = 34) 
     for (seed in seed.vec) {
          load(file = paste("data-sanity/crt", n1, "_seed", seed, "_", st, ".RData", sep = ""))
          mle_tmp_dat = cbind(mle_tmp_dat, est.mat)
     }
     mle.mat = mle_tmp_dat[,-1]
     
     alpha.idx = c(1,c(6 * (0:4) + 5))
     alpha.sd.idx = c(3, c(6 * (0:4) + 8))
     
     bias = apply(mle.mat[alpha.idx, ], 1, mean) - rep(gamma.true, length(alpha.idx))
     se = apply(mle.mat[alpha.idx, ], 1, sd) / sqrt(1000)
     
     bias.sd = paste(round(bias*100, 3), "(", round(se*100, 3), ")", sep = "")
     
     mcsd = apply(mle.mat[alpha.idx, ], 1, sd)
     
     sd = apply(mle.mat[alpha.sd.idx, ], 1, mean)
     
     accuracy = round(sd / mcsd, 3)
     
     lowerB = mle.mat[alpha.idx,] - 1.96 * mle.mat[alpha.sd.idx,]
     upperB = mle.mat[alpha.idx,] + 1.96 * mle.mat[alpha.sd.idx,]
     
     #cr = round(apply(((rep(gamma.true, 5) > lowerB) * (rep(gamma.true, 5) < upperB)), 1, sum) / 1000,3)
     cr = round(apply(((rep(gamma.true, length(alpha.idx)) > lowerB) * (rep(gamma.true, length(alpha.idx)) < upperB)), 1, mean),3)
     
     rst.tmp = matrix(c(bias.sd, accuracy, cr), ncol = 3)
     
     rownames(rst.tmp) = c(paste("cohort", n1, sep = ""), paste("cohort", n1, "_casecrt", n2.vec, sep = ""))
     
     colnames(rst.tmp) = c("bias_sd", "acc", "coverage")
     final = rbind(final, rst.tmp)
}
final[-1,]
