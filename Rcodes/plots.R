rm(list = ls())
library(brm)
library(tidyverse)
library(latex2exp)

beta = c(-7, 0.5)
gamma = c(0,-1)
eta = c(0.2,-1)

v = seq(from = -2, to = 2, by = 0.01)
va = cbind(rep(1, length(v)),v)

# pr(z=1|v) --------------
pi = exp(va %*% eta) / (1 + exp(va %*% eta))

# log(rr) ----------------
logrr = va %*% gamma

# log(op) ---------------
logop = va %*% beta

# pr(y=1|v) -----------------
P.mat = t(mapply(getProbScalarRR, va %*% gamma, va %*% beta))
q = apply(P.mat * cbind(1-pi, pi), 1, sum)


dat = data.frame(v, pi, q, logrr, logop)

# plots -------------------
theme_set(theme_bw() + 
               theme(axis.text.x = element_text(color = "grey20", size = 18, face = "plain"),
                     axis.text.y = element_text(color = "grey20", size = 18, face = "plain"),
                     axis.title.x = element_text(color = "grey20", size =20, face = "bold"),
                     axis.title.y = element_text(color = "grey20", size = 20, face = "bold"),
                     legend.text = element_text(size = 18),
                     legend.title = element_text(size = 20, face = "bold"),
                     strip.text.x = element_text(size=18, face="bold")))

dat %>% 
     ggplot(aes(x = v, y = pi)) +
     geom_line(size = 1) + 
     ylab(TeX('$\\Pr(Z=1|v_1)$')) + 
     xlab(TeX('$\\v_1$')) + 
     scale_y_continuous(breaks = seq(0,1,0.2))
ggsave("figs/pi.png", width = 6, height = 4)

dat %>% 
     ggplot(aes(x = v, y = q)) +
     geom_line(size = 1) + 
     ylab(TeX('$\\Pr(Y=1|v_1)$')) + 
     xlab(TeX('$\\v_1$')) 
ggsave("figs/pry1v.png", width = 6, height = 4)

dat %>% 
     ggplot(aes(x = v, y = logrr)) +
     geom_line(size = 1) + 
     ylab(TeX("$\\log RR(v_1)$"))+
     xlab(TeX('$\\v_1$')) + 
ggsave("figs/logrr.png", width = 6, height = 4)

dat %>% 
     ggplot(aes(x = v, y = logop)) +
     geom_line(size = 1) + 
     ylab(TeX("$\\log OP(v_1)$"))+
     xlab(TeX('$\\v_1$')) + 
ggsave("figs/logop.png", width = 6, height = 4)
                                                                                                                                        
# histgram for v|Y=1 --------------
source("Rcodes/0c_simulation.R")
casectrl.dat = simCasecrt(n = 100000, pr = 1/2, gamma, beta, eta)

casectrl.dat %>%
     filter(y == 1) %>%
     ggplot(aes(x = v2))  + 
     geom_histogram(aes(y=..density..), colour="black", fill="white")+
     geom_density(alpha=.2, fill="#FF6666") +
     xlab(TeX('$\\v_1$')) + ylab(TeX('$\\Pr(v_1|Y=1)$'))
ggsave("figs/density_vy1.png", width = 6, height = 4)

casectrl.dat %>%
     filter(y == 0) %>%
     ggplot(aes(x = v2))  + 
     geom_histogram(aes(y=..density..), colour="black", fill="white")+
     geom_density(alpha=.2, fill="#FF6666") +
     xlab(TeX('$\\v_1$')) + ylab(TeX('$\\Pr(v_1|Y=0)$'))
ggsave("figs/density_vy0.png", width = 6, height = 4)
     