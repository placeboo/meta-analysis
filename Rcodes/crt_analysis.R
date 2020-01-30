rm(list = ls())

n.vec = c(500, 1000, 5000)

st = "01292020"

seed.vec = 1:10

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)

rst.dat = matrix(NA, ncol = 6)
for (n in n.vec) {
     
     mle_tmp_dat = matrix(NA, nrow = 8)
     for (seed in seed.vec) {
          load(file = paste("crt/crt", n,"_seed", seed, "_", st, ".RData", sep = ""))
          mle_tmp_dat = cbind(mle_tmp_dat, est.mat)
     }
     mle.mat = mle_tmp_dat[,-1]
     
     bias = apply(mle.mat[c(1, 2), ], 1, mean) - gamma.true
     se = apply(mle.mat[c(1, 2), ], 1, sd) / sqrt(1000)
     
     bias.sd = paste(round(bias,3), "(", round(se,3), ")", sep = "")
     
     mcsd = apply(mle.mat[c(1, 2), ], 1, sd)
     
     sd = apply(mle.mat[c(5,6), ], 1, mean)
     
     accuracy = round(sd / mcsd, 3)
     
     lowerB = mle.mat[c(1, 2),] - 1.96 * mle.mat[c(5,6),]
     upperB = mle.mat[c(1, 2),] + 1.96 * mle.mat[c(5,6),]
     
     cr = round(apply(((gamma.true > lowerB) * (gamma.true < upperB)), 1, mean)  ,3)
     
     
     rst.tmp = cbind(matrix(bias.sd, ncol = 2, byrow = TRUE),
                     matrix(accuracy, ncol = 2, byrow = TRUE),
                     matrix(cr, ncol = 2, byrow = TRUE))
     
     rst.dat = rbind(rst.dat, rst.tmp) 
}
rst.dat[-1,]
