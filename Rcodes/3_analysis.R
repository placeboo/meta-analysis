rm(list = ls())

n1.vec = c(500, 1000, 5000)
n2.vec = n1.vec
st = "01282020"

seed.vec = 1:10

beta.true = c(-7, 0.5)
gamma.true = c(0,-1)
eta.true = c(0.2,-1)

bias.mat = matrix(NA,  ncol = 2)
accuracy.mat = matrix(NA,  ncol = 2)
cr.mat = matrix(NA,  ncol = 2)
na.count = rep(NA)

final = matrix(NA, ncol = 6)

for (n1 in n1.vec) {
     rst.dat = matrix(NA, ncol = 6)
     for (n2 in 500) {
          mle_tmp_dat = matrix(NA, nrow = 40)
          for (seed in seed.vec) {
               load(file = paste("data/crt", n1, "_casecrt", n2, "_seed", seed, "_", st, ".RData", sep = ""))
               mle_tmp_dat = cbind(mle_tmp_dat, est.mat)
          }
          mle.mat = mle_tmp_dat[,-1]
          
          bias = apply(mle.mat[c(1, 2, 9, 10, 17, 18, 25, 26, 33, 34), ], 1, mean) - rep(gamma.true, 5)
          se = apply(mle.mat[c(1, 2, 9, 10, 17, 18, 25, 26, 33, 34), ], 1, sd) / sqrt(1000)
          
          bias.sd = paste(round(bias,3), "(", round(se,3), ")", sep = "")
          
          mcsd = apply(mle.mat[c(1, 2, 9, 10, 17, 18, 25, 26, 33, 34), ], 1, sd)
          
          sd = apply(mle.mat[c(5,6,13,14,21,22,29,30,37,38), ], 1, mean)
          
          accuracy = round(sd / mcsd, 3)
          
          lowerB = mle.mat[c(1, 2, 9, 10, 17, 18, 25, 26, 33, 34),] - 1.96 * mle.mat[c(5,6,13,14,21,22,29,30,37,38),]
          upperB = mle.mat[c(1, 2, 9, 10, 17, 18, 25, 26, 33, 34),] + 1.96 * mle.mat[c(5,6,13,14,21,22,29,30,37,38),]
          
          cr = round(apply(((rep(gamma.true, 5) > lowerB) * (rep(gamma.true, 5) < upperB)), 1, sum) / 1000,3)
          
           
          rst.tmp = cbind(matrix(bias.sd, ncol = 2, byrow = TRUE),
                matrix(accuracy, ncol = 2, byrow = TRUE),
                matrix(cr, ncol = 2, byrow = TRUE))
          
          rst.dat = rbind(rst.dat, rst.tmp)
     }
     rst.dat = rst.dat[-1, ]
     final = rbind(final, rst.dat)
}
