rm(list = ls())

n1.vec = c(500, 1000, 2000, 5000, 10000)
n2.vec = c(200, 500, 1000, 5000)
st = "02042020"

seed.vec = 1:10

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)

# beta.true = c(-1, 0.5)
# gamma.true = c(0,1)
# eta.true = c(1,-1)

final = matrix(NA, ncol = 6)

for (n1 in n1.vec) {
     rst.dat = matrix(NA, ncol = 6)
   
     for (n2 in n2.vec) {
          mle_tmp_dat = matrix(NA, nrow = 28)
          for (seed in seed.vec) {
               load(file = paste("data4/crt", n1, "_casecrt", n2, "_seed", seed, "_", st, ".RData", sep = ""))
               mle_tmp_dat = cbind(mle_tmp_dat, est.mat)
          }
          mle.mat = mle_tmp_dat[,-1]
          
          bias = apply(mle.mat[c(1, 2, 9, 10, 17, 18), ], 1, mean) - rep(gamma.true, 3)
          se = apply(mle.mat[c(1, 2, 9, 10, 17, 18), ], 1, sd) / sqrt(1000)
          
          bias.sd = paste(round(bias*100, 3), "(", round(se*100, 3), ")", sep = "")
          
          mcsd = apply(mle.mat[c(1, 2, 9, 10, 17, 18), ], 1, sd)
          
          sd = apply(mle.mat[c(5, 6, 13, 14, 23, 24), ], 1, mean)
          
          accuracy = round(sd / mcsd, 3)
          
          lowerB = mle.mat[c(1, 2, 9, 10, 17, 18),] - 1.96 * mle.mat[c(5, 6, 13, 14, 23, 24),]
          upperB = mle.mat[c(1, 2, 9, 10, 17, 18),] + 1.96 * mle.mat[c(5, 6, 13, 14, 23, 24),]
          
          #cr = round(apply(((rep(gamma.true, 5) > lowerB) * (rep(gamma.true, 5) < upperB)), 1, sum) / 1000,3)
          cr = round(apply(((rep(gamma.true, 3) > lowerB) * (rep(gamma.true, 3) < upperB)), 1, mean),3)
           
          rst.tmp = cbind(matrix(bias.sd, ncol = 2, byrow = TRUE),
                matrix(accuracy, ncol = 2, byrow = TRUE),
                matrix(cr, ncol = 2, byrow = TRUE))
          
          rst.dat = rbind(rst.dat, rst.tmp)
     }
     rst.dat = rst.dat[-1, ]
     final = rbind(final, rst.dat)
}
final[-1,]
