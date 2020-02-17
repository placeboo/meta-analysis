rm(list = ls())

n1.vec = c(500, 1000, 2000, 4000, 8000)
n2.vec = c(200, 500, 1000, 2000, 4000)
st = "02082020"

seed.vec = 1:50

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)


final = matrix(NA, ncol = 6)

for (n1 in n1.vec) {
     mle_tmp_dat = matrix(NA, nrow = 108)
     for (seed in seed.vec) {
          filename = paste("data6/crt", n1, "_seed", seed, "_", st, ".RData", sep = "")
          if(file.exists(filename)) {
               load(filename)
          } else {
               next
          }
          mle_tmp_dat = cbind(mle_tmp_dat, est.mat)
     }
     mle.mat = mle_tmp_dat[,-1]
     
     alpha.idx = sort(c(c(1, 2),
                        c(20 * (0:4) + 9), 
                        c(20 * (0:4) + 10), 
                        c(17 + 20 * (0:4)), 
                        c(18 + 20 * (0:4))))
     alpha.sd.idx = sort(c(c(5,6), 
                           c(13 + 20 * (0:4)),
                           c(14 + 20 * (0:4)),
                           c(23 + 20 * (0:4)),
                           c(24 + 20 * (0:4))))
     bias = apply(mle.mat[alpha.idx, ], 1, mean) - rep(gamma.true, length(alpha.idx)/2)
     se = apply(mle.mat[alpha.idx, ], 1, sd) / sqrt(ncol(mle.mat))
     
     bias.sd = paste(round(bias*100, 3), "(", round(se*100, 3), ")", sep = "")
     
     mcsd = apply(mle.mat[alpha.idx, ], 1, sd)
     
     sd = apply(mle.mat[alpha.sd.idx, ], 1, mean)
     
     accuracy = round(sd / mcsd, 3)
     
     lowerB = mle.mat[alpha.idx,] - 1.96 * mle.mat[alpha.sd.idx,]
     upperB = mle.mat[alpha.idx,] + 1.96 * mle.mat[alpha.sd.idx,]
     
     #cr = round(apply(((rep(gamma.true, 5) > lowerB) * (rep(gamma.true, 5) < upperB)), 1, sum) / 1000,3)
     cr = round(apply(((rep(gamma.true, length(alpha.idx)/2) > lowerB) * (rep(gamma.true, length(alpha.idx)/2) < upperB)), 1, mean),3)
     
     rst.tmp = cbind(matrix(bias.sd, ncol = 2, byrow = TRUE),
                     matrix(accuracy, ncol = 2, byrow = TRUE),
                     matrix(cr, ncol = 2, byrow = TRUE))
     
     rowname = rep(NA, nrow(rst.tmp))
     rowname[1] = paste("cohort", n1, sep = "")
     rowname[seq(from = 2, to = nrow(rst.tmp), by = 2)] = paste("casecrt", n2.vec, sep = "")
     rowname[seq(from = 3, to = nrow(rst.tmp), by = 2)] = paste("cohort", n1, "_casecrt", n2.vec, sep = "")
     
     rownames(rst.tmp) = rowname
     colnames(rst.tmp) = c(paste("bias_sd", c(1,2), sep = ""), paste("acc", c(1,2), sep = ""), paste("coverage", c(1,2), sep = ""))
     final = rbind(final, rst.dat)
}
final[-1,]
