rm(list = ls())

# n = 2000
# n2.vec = seq(from = 0, to = 1600, by=200)

# n = 3000
# n2.vec = seq(from = 0, to = 2400, by = 300)
# n = 4000
# n2.vec = seq(from = 0, to=3600, by = 400)

# n = 8000
# n2.vec = seq(from = 0, to = 7200, by = 800)

n = 10000
n2.vec = seq(from = 0, to = 9000, by = 1000)

# n = 16000
# n2.vec = seq(from = 0, to = 14400, by = 1600)
st = "02152020"

seed.vec = 1:25

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)


final = matrix(NA, ncol = 6)

for (n2 in n2.vec) {
     n1 = n - n2
     
     if (n2 == 0) {
          mle_tmp_dat = matrix(NA, nrow = 8)
     } else {
          mle_tmp_dat = matrix(NA, nrow = 12)
     }
     
     for (seed in seed.vec) {
          # filename = paste("data6/crt", n1, "_seed", seed, "_", st, ".RData", sep = "")
          # if(file.exists(filename)) {
          #      load(filename)
          # } else {
          #      next
          # }
          
          filename = paste("data7/crt", n1, "_casecrt", n2,"_seed", seed, "_", st, ".RData", sep = "")
          load(filename)
          mle_tmp_dat = cbind(mle_tmp_dat, est.mat)
     }
     mle.mat = mle_tmp_dat[,-1]
     
     alpha.idx = c(1,2)
     if (n2 == 0) {
          alpha.sd.idx = c(5,6)
     } else{
          alpha.sd.idx = c(7,8)
     }
     
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
     rownames(rst.tmp) = paste("cohort", n1, "_casecrt", n2, sep = "")
     colnames(rst.tmp) = c(paste("bias_sd", c(1,2), sep = ""), paste("acc", c(1,2), sep = ""), paste("coverage", c(1,2), sep = ""))
     final = rbind(final, rst.tmp)
}
final[-1,]
