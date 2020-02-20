rm(list = ls())

seed.vec = 1:5
n.vec = c(500,1000,2000, 4000, 8000)
st = "02182020"

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)

rst.dat = matrix(NA, ncol = 6)

for (n in n.vec) {
     mle_tmp_dat = matrix(NA, nrow = 12)
     for (seed in seed.vec) {
          load(file = paste("data-case-control/casecrt", n, "_seed", seed, "_", st, ".RData", sep = ""))
          mle_tmp_dat = cbind(mle_tmp_dat, est.mat)
     }
     mle.mat = mle_tmp_dat[,-1]
     
     alpha.idx = c(1,2)
     alpha.sd.idx = c(7,8)
     bias = apply(mle.mat[alpha.idx, ], 1, mean) - rep(gamma.true, length(alpha.idx)/2)
     se = apply(mle.mat[alpha.idx, ], 1, sd) / sqrt(ncol(mle.mat))
     
     bias.sd = paste(round(bias*100, 3), "(", round(se*100, 3), ")", sep = "")
     
     mcsd = apply(mle.mat[alpha.idx, ], 1, sd)
     
     sd = apply(mle.mat[alpha.sd.idx, ], 1, mean)
     
     accuracy = round(sd / mcsd, 3)
     
     lowerB = mle.mat[alpha.idx,] - 1.96 * mle.mat[alpha.sd.idx,]
     upperB = mle.mat[alpha.idx,] + 1.96 * mle.mat[alpha.sd.idx,]
     
     cr = round(apply(((rep(gamma.true, length(alpha.idx)/2) > lowerB) * (rep(gamma.true, length(alpha.idx)/2) < upperB)), 1, mean),3)
     
     rst.tmp = cbind(matrix(bias.sd, ncol = 2, byrow = TRUE),
                     matrix(accuracy, ncol = 2, byrow = TRUE),
                     matrix(cr, ncol = 2, byrow = TRUE))
     rownames(rst.tmp) = paste("casecrt", n, sep = "")
     colnames(rst.tmp) = c(paste("bias_sd", c(1,2), sep = ""), paste("acc", c(1,2), sep = ""), paste("coverage", c(1,2), sep = ""))
     
     rst.dat = rbind(rst.dat, rst.tmp) 
}
rst.dat[-1,]
