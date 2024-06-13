#vaje2
library(matlib)
library(dplyr)
podatki <- read.csv("data2.csv")

X <- as.matrix(podatki[2])
Z <- as.matrix(podatki[-c(1,2)])

norma2 <- function(v){
  sqrt(sum(v^2))
}

grahm_schmidt <- function(Z, eps) {
  m <- nrow(Z)
  d <- ncol(Z)
  S <- matrix(0, nrow = m, ncol = d)
  P <- matrix(0, nrow = d, ncol = d)
  for (j in 1:d) {
    v <- Z[ , j, drop = FALSE]
    if (j > 1) {
      for(i in 1:(j-1)) {
        P[i, j] <- t(S[,i,drop = FALSE]) %*% Z[ , j, drop = FALSE]
        v <- v - P[i, j] * S[ ,i]
      }
    }
    if (norma2(v) < eps){
      S[ ,j] <- 0
      P[j,j] <- 1
    }
    else{
      P[j, j] = norma2(v)
      S[, j] = v/ P[j, j]
    }
  }
  
  return(list(S=S, P=P))
  
}

eps <- 10^-10

razcep <- grahm_schmidt(Z, eps)
S <- razcep$S
P <- razcep$P


Z_nov <- S %*% P 


norme_S <- apply(S, 2, norma2)
J <- diag(norme_S)

ztz <- t(P) %*% J %*% P
ztz_inv <- inv(P) %*% J %*% inv(t(P))

beta_hat <- ztz_inv %*% t(Z) %*% X


# 3 vaje - linearna regresija
primer2 <- read.csv("primer2.csv")
sigma <- 1
ei <- rnorm(300, 0, sigma^2)

Z.2 <- primer2[-c(1,2)]
Z.2 <- as.matrix(cbind(z0 = 1, Z.2))
X.2 <- as.matrix(primer2[2])

# prvi način - GS algoritem

razcep.2 <- grahm_schmidt(Z.2, eps)
S.2 <- razcep.2$S
P.2 <- razcep.2$P

norme_S.2 <- apply(S.2, 2, norma2)
J.2 <- diag(norme_S.2)

ztz.2 <- t(P.2) %*% J.2 %*% P.2
ztz_inv.2 <- solve(P.2) %*% J.2 %*% solve(t(P.2))

beta_hat.2.posploseni <- ztz_inv.2 %*% t(Z.2) %*% X.2



# drugi način - navaden inverz

beta_hat.2.navadni <- solve(t(Z.2) %*% Z.2) %*% t(Z.2) %*% X.2


# tretji način - z lm

lm.model <- lm(X.2 ~ Z.2) 
beta_hat.2.lm <- summary(lm.model)$coeff
summary(lm.model)

# 3. naloga - ANOVA z enojno klasifikacijo

m <- 4
mu <- 2
alphe <- c(0.2, -0.5, 1, 2)
sigme <- c(0.05, 1, 0.03, 0.08)
nji <- c(25, 30, 25, 27)
n <- sum(nji)
sigma <- rep(0, n)

Z.3 <- matrix(0, nrow = n, ncol = m + 1)
Z.3[,1] <- 1
zacetek <- 1
for (j in 1:m){
  Z.3[zacetek:(zacetek+nji[j]-1), j+1] <- 1
  sigma[zacetek:(zacetek+nji[j]-1)] <- sigme[j]
  zacetek <- zacetek + nji[j]
}

epsi <- rnorm(n, 0, sigma^2)
X.3 <- Z.3 %*% alpha_mu + epsi


#GS
razcep.3 <- grahm_schimdt(Z.3, eps)
S.3 <- razcep.3$S
P.3 <- razcep.3$P

norme_S.3 <- apply(S.3, 2, norma2)
J.3 <- diag(norme_S.3)

ztz.3 <- t(P.3) %*% J.3 %*% P.3
ztz_inv.3 <- solve(P.3) %*% J.3 %*% solve(t(P.3))

parametri <- ztz_inv.3 %*% t(Z.3) %*% X.3
