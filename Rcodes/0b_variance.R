var.1d = function(gamma, beta, eta, casectrl.dat, cohort.dat) {
     
     cohort.va = as.matrix(cohort.dat$va)
     cohort.vb = cohort.va
     cohort.vc = cohort.va
     cohort.Pmat = t(mapply(getProbScalarRR, cohort.va %*% gamma, cohort.vb %*% beta))
     cohort.pi = exp(cohort.vc %*% eta) / (exp(cohort.vc %*% eta) + 1)    
     cohort.p.vec = rep(NA, nrow(cohort.va))
     for (j in 1:nrow(cohort.dat)){
          cohort.p.vec[j] = cohort.Pmat[j, cohort.dat$z[j]+1]
     }
     
     pl1.pgamma = (cohort.dat$y - cohort.p.vec) / (1 - cohort.p.vec) * 
          (cohort.dat$z - (1-cohort.Pmat[,1])/(1-cohort.Pmat[,1] + 1-cohort.Pmat[,2])) * cohort.va
     
     pl1.pbeta =  (cohort.dat$y - cohort.p.vec) / (1 - cohort.p.vec) * ((1-cohort.Pmat[,1]) * (1 - cohort.Pmat[,2])) / ((1-cohort.Pmat[,1]) + (1 - cohort.Pmat[,2])) * cohort.vb
     
     pl1.peta = 0 
     
     # l2 deviatives -------------
     pl2.pgamma = 0
     pl2.pbeta = 0
     pl2.peta = (cohort.dat$z - cohort.pi) * cohort.vc 
     
     # l3 deviatives ------------
     casectrl.va = as.matrix(casectrl.dat$va)
     casectrl.vb = casectrl.va
     casectrl.vc = casectrl.va
     casectrl.Pmat = t(mapply(getProbScalarRR, casectrl.va %*% gamma, casectrl.vb %*% beta))
     casectrl.pi = exp(casectrl.vc %*% eta) / (exp(casectrl.vc %*% eta) + 1)
     
     casectrl.p.vec = rep(NA, nrow(casectrl.dat))
     for (j in 1:nrow(casectrl.dat)){
          casectrl.p.vec[j] = casectrl.Pmat[j, casectrl.dat$z[j]+1]
     }
     qA = (casectrl.dat$y *  casectrl.Pmat[,2] + (1 - casectrl.dat$y) * (1 -  casectrl.Pmat[,2])) * casectrl.pi
     qB = qA + (casectrl.dat$y *  casectrl.Pmat[,1] + (1 - casectrl.dat$y) * (1 -  casectrl.Pmat[,1])) * (1 - casectrl.pi)
     casectrl.q = qA/qB 
     
     tmp = (casectrl.dat$z - casectrl.q) / (casectrl.q * (1 -casectrl.q))
     
     pp1.pgamma = casectrl.Pmat[,2] * (1 - casectrl.Pmat[,2])/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2])) * casectrl.va
     pp0.pgamma = casectrl.Pmat[,1] * ( - (1 - casectrl.Pmat[,1])/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.va
     
     pp1.pbeta = casectrl.Pmat[,2] * (((1 - casectrl.Pmat[,1]) * (1 - casectrl.Pmat[,2]))/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.vb
     pp0.pbeta = casectrl.Pmat[,1] * (((1 - casectrl.Pmat[,1]) * (1 - casectrl.Pmat[,2]))/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.vb
     
     pA.pgamma = (2 * casectrl.dat$y - 1)  * pp1.pgamma * casectrl.pi
     pA.pbeta = (2 * casectrl.dat$y - 1)  * pp1.pbeta * casectrl.pi
     pA.peta = (casectrl.dat$y * casectrl.Pmat[,2] + (1 - casectrl.dat$y) * (1 - casectrl.Pmat[,2])) * casectrl.pi * (1 - casectrl.pi) * casectrl.vc
     
     pB.pgamma = pA.pgamma + (2 * casectrl.dat$y - 1) * pp0.pgamma * (1-casectrl.pi)
     pB.pbeta = pA.pbeta + (2 * casectrl.dat$y - 1) * pp0.pbeta * (1-casectrl.pi)
     pB.peta = pA.peta -  (casectrl.dat$y * casectrl.Pmat[,1] + (1 - casectrl.dat$y) * (1 - casectrl.Pmat[,1])) *  casectrl.pi * (1 - casectrl.pi) * casectrl.vc
     
     pq.pgamma = (qB * pA.pgamma - qA * pB.pgamma) * 1/qB * 1/qB
     pq.pbeta = (qB * pA.pbeta - qA * pB.pbeta) * 1/qB * 1/qB
     pq.peta = (qB * pA.peta- qA * pB.peta) * 1/qB * 1/qB
     
     pl3.pgamma = tmp * pq.pgamma
     pl3.pbeta = tmp * pq.pbeta
     pl3.peta = tmp * pq.peta
     
     # combina
     dl.dgamma = c(pl1.pgamma, pl3.pgamma)
     dl.dbeta = c(pl1.pbeta, pl3.pbeta)
     dl.deta = c(pl2.peta, pl3.peta)
     
     dl = cbind(dl.dgamma, dl.dbeta, dl.deta)
     hessian.mat = t(dl) %*% dl
     var.mat = solve(hessian.mat)
     return(var.mat)
}


var.meta = function(gamma, beta, eta,
                    case.control = list(y, z, va, vb, vc),
                    cohort = list(y, z, va, vb, vc)) {
     
     cohort.va = cohort$va
     cohort.vb = cohort$vb
     cohort.vc = cohort$vc
     
     cohort.Pmat = t(mapply(getProbScalarRR, cohort.va %*% gamma, cohort.vb %*% beta))
     
     cohort.pi = c(exp(cohort.vc %*% eta) / (exp(cohort.vc %*% eta) + 1))
     cohort.p.vec = rep(NA, nrow(cohort.va))
     for (j in 1:nrow(cohort.dat)){
          cohort.p.vec[j] = cohort.Pmat[j, cohort.dat$z[j]+1]
     }
     
     pl1.pgamma = (cohort.dat$y - cohort.p.vec) / (1 - cohort.p.vec) * 
          (cohort.dat$z - (1-cohort.Pmat[,1])/(1-cohort.Pmat[,1] + 1-cohort.Pmat[,2])) * cohort.va
     
     pl1.pbeta =  (cohort.dat$y - cohort.p.vec) / (1 - cohort.p.vec) * ((1-cohort.Pmat[,1]) * (1 - cohort.Pmat[,2])) / ((1-cohort.Pmat[,1]) + (1 - cohort.Pmat[,2])) * cohort.vb
     
     pl1.peta = 0 
     
     # l2 deviatives -------------
     pl2.pgamma = 0
     pl2.pbeta = 0
     pl2.peta = (cohort.dat$z - cohort.pi) * cohort.vc 
     
     # l3 deviatives ------------
     casectrl.va = case.control$va
     casectrl.vb = case.control$vb
     casectrl.vc = case.control$vc
     casectrl.Pmat = t(mapply(getProbScalarRR, casectrl.va %*% gamma, casectrl.vb %*% beta))
     casectrl.pi = c(exp(casectrl.vc %*% eta) / (exp(casectrl.vc %*% eta) + 1))
     
     casectrl.p.vec = rep(NA, nrow(casectrl.dat))
     for (j in 1:nrow(casectrl.dat)){
          casectrl.p.vec[j] = casectrl.Pmat[j, casectrl.dat$z[j]+1]
     }
     qA = c((casectrl.dat$y *  casectrl.Pmat[,2] + (1 - casectrl.dat$y) * (1 -  casectrl.Pmat[,2])) * casectrl.pi)
     qB = qA + c((casectrl.dat$y *  casectrl.Pmat[,1] + (1 - casectrl.dat$y) * (1 -  casectrl.Pmat[,1])) * (1 - casectrl.pi))
     casectrl.q = qA/qB 
     
     tmp = (casectrl.dat$z - casectrl.q) / (casectrl.q * (1 -casectrl.q))
     
     pp1.pgamma = casectrl.Pmat[,2] * (1 - casectrl.Pmat[,2])/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2])) * casectrl.va
     pp0.pgamma = casectrl.Pmat[,1] * ( - (1 - casectrl.Pmat[,1])/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.va
     
     pp1.pbeta = casectrl.Pmat[,2] * (((1 - casectrl.Pmat[,1]) * (1 - casectrl.Pmat[,2]))/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.vb
     pp0.pbeta = casectrl.Pmat[,1] * (((1 - casectrl.Pmat[,1]) * (1 - casectrl.Pmat[,2]))/((1 - casectrl.Pmat[,1]) + (1 - casectrl.Pmat[,2]))) * casectrl.vb
     
     pA.pgamma = (2 * casectrl.dat$y - 1)  * casectrl.pi * pp1.pgamma 
     pA.pbeta = (2 * casectrl.dat$y - 1) * casectrl.pi * pp1.pbeta 
     pA.peta = (casectrl.dat$y * casectrl.Pmat[,2] + (1 - casectrl.dat$y) * (1 - casectrl.Pmat[,2])) * casectrl.pi * (1 - casectrl.pi) * casectrl.vc
     
     pB.pgamma = pA.pgamma + (2 * casectrl.dat$y - 1) * pp0.pgamma * (1-casectrl.pi)
     pB.pbeta = pA.pbeta + (2 * casectrl.dat$y - 1) * pp0.pbeta * (1-casectrl.pi)
     pB.peta = pA.peta -  (casectrl.dat$y * casectrl.Pmat[,1] + (1 - casectrl.dat$y) * (1 - casectrl.Pmat[,1])) *  casectrl.pi * (1 - casectrl.pi) * casectrl.vc
     
     pq.pgamma = (qB * pA.pgamma - qA * pB.pgamma) * 1/qB * 1/qB
     pq.pbeta = (qB * pA.pbeta - qA * pB.pbeta) * 1/qB * 1/qB
     pq.peta = (qB * pA.peta- qA * pB.peta) * 1/qB * 1/qB
     
     pl3.pgamma = tmp * pq.pgamma
     pl3.pbeta = tmp * pq.pbeta
     pl3.peta = tmp * pq.peta
     
     # combina
     dl.dgamma = rbind(pl1.pgamma, pl3.pgamma)
     dl.dbeta = rbind(pl1.pbeta, pl3.pbeta)
     dl.deta = rbind(pl2.peta, pl3.peta)
     
     dl = cbind(dl.dgamma, dl.dbeta, dl.deta)
     hessian.mat = t(dl) %*% dl
     var.mat = solve(hessian.mat)
     return(var.mat)
}
