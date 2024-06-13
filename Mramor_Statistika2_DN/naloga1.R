# 1. naloga

# uvozimo podatke
podatki1 <- read.csv("podatki_1.txt", header=FALSE)


# a) Poiščite kompletno zadostno statistiko

T_1 <- sum(log(podatki1$V1))
T_2 <- sum(podatki1$V1)

##############################################################################

# b) Cenilka po metodi največjega verjetja
# zapišemo g kot funkcijo in rešujemo z Newton-Raphsonom numerično - iščemo ničlo funkcije g

g <- function(x) digamma(x) - log(x) - mean(log(podatki1$V1)) + log(mean(podatki1$V1))

g_odvod <- function(x) trigamma(x) - 1/x

a_cnv <- 0.1
for(i in 1:10000){
  a_cnv <- a_cnv - g(a_cnv)/g_odvod(a_cnv)
}

# preverimo če smo res našli ničlo
g(a_cnv)

# rezultat
a_cnv 

# alternativno, brez uporabe for zanke
a_cnv_alt <- uniroot(g, c(0.00001, a_cnv + 4))$root
a_cnv_alt

# poračunamo še b
b_cnv = a_cnv/mean(podatki1$V1)
b_cnv

##############################################################################

# c) Cenilka po metodi momentov

mu1 <- mean(podatki1$V1)
mu2 <- mean(podatki1$V1^2)
a_mm <- mu1^2/(mu2 - mu1^2)
a_mm
b_mm <- mu1/(mu2 - mu1^2)
b_mm

#############################################################################

# d) Metoda delta - poskus implementacije v R

library(msm)

podatki1$x2 <- podatki1$V1^2
colnames(podatki1) <- c("x1", "x2")

cov_est <- cov(podatki1)

# metoda delta vrne aproksimacijo standardne napake transformacije g(X) za 
# slučanjo spremenljivko X=(x_1,x_2,...), če ji podamo prvih m-momentov

# za (m razsežno funkcijo g) in variančno-kovariančno matriko X
deltamethod(list(~ x1^2/(x2-x1^2), ~ x1/(x2-x1^2)), c(mu1, mu2), cov_est)

################DODANO######################################################

Jacobi <- matrix(c(2*mu1*mu2/(mu2-mu1^2)^2,
                   (mu2+mu1^2)/(mu2-mu1^2)^2,
                   - mu1^2/(mu2-mu1^2)^2,
                   - mu1/(mu2-mu1^2)^2),
                 ncol=2, nrow=2)
JacobiT <- t(Jacobi)
sigma <- matrix(c(mu2-mu1^2,
                  2*mu2*(mu2-mu1^2)/mu1,
                  2*mu2*(mu2-mu1^2)/mu1,
                  mu2*(6*mu2-2*mu1^2)*(mu2-mu1^2)/mu1^2),
                ncol=2, nrow=2)

SigmaD <- Jacobi %*% sigma %*% JacobiT

# prepisana iz računa na roke
SigmaD_2 <- matrix(c(2*mu1^2*mu2/(mu1^4+mu2^2-2*mu1^2*mu2),
                     2*mu1*mu2/(mu1^4+mu2^2-2*mu1^2*mu2),
                     2*mu1*mu2/(mu1^4+mu2^2-2*mu1^2*mu2),
                     (-mu1^2+3*mu2)/(mu1^4+mu2^2-2*mu1^2*mu2)),
                   ncol=2, nrow=2)
# res pride enako, tako da izgleda vredu

##########################################################################

# e) Aproksimativno območje zaupanja stopnje zaupanja 0.95

n=nrow(podatki1)
alpha = 0.95
d = 2
fisher = qf(alpha,d,n-d)
ksi2 = d*(n-1)/(n-d)*fisher

# izračunamo inverz matrike SigmaD/n in jo diagonaliziramo
inv_sigma = solve(SigmaD/n)
lambda = diag(eigen(inv_sigma)$values)
Q = eigen(inv_sigma)$vectors

# območje zaupanja je elipsoid z polosmi dolžine:
dolzina1 <- sqrt(ksi2/lambda[1,1])
dolzina1
dolzina2 <- sqrt(ksi2/lambda[2,2])
dolzina2
