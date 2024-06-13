# 5. naloga

podatki5 <- read.delim("podatki_5.txt")

# iz vaj imamo funkciji za normo in GS algoritem

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

##############################################################################

# a) Predpostavimo linerani model:
# LDL = Beta_0 + Beta_1 * TCH + Beta_2 * HDL + Beta_3 * TRI + epsilon

# (i) Ocenite Beta0, Beta1, Beta2, Beta3 po MNK

# definiramo Z, X in epsilon = ei, z sigma = 1
sigma <- 1
ei <- rnorm(300, 0, sigma^2)

Z <- podatki5[c(1:3)]
Z <- as.matrix(cbind(Z0 = 1, Z))
X <- as.matrix(podatki5[4])

# prvi način - GS algoritem

#nastavimo mejo
eps <- 10^-10 
razcep <- grahm_schmidt(Z, eps)
S <- razcep$S
P <- razcep$P

norme_S <- apply(S, 2, norma2)
J <- diag(norme_S)

ztz <- t(P) %*% J %*% P
ztz_inv <- solve(P) %*% J %*% solve(t(P))

beta_hat_posploseni <- ztz_inv %*% t(Z) %*% X


# drugi način - navaden inverz

beta_hat_navadni <- solve(t(Z) %*% Z) %*% t(Z) %*% X


# tretji način - z lm

lm.model <- lm(data=podatki5, LDL ~ TCH + HDL + TRI) 
summary(lm.model)
beta_hat_lm <- summary(lm.model)$coeff
#############################################################################

# (ii) Preizkusite domnevo Beta1=Beta2=Beta3=0 na standarden način

# ne gre z razmerjem verjetij kot na vajah zato postopamo s Fisherjevim F testom.

# poračunamo VKR za unrestricted model
VKR_u <- sum(resid(lm.model)^2)

# izračunamo še restricted model in poračunamo VKR
restricted.lm.model <- lm(data=podatki5, LDL ~ 1)
summary(restricted.lm.model)
VKR_r <- sum(resid(restricted.lm.model)^2)

# nastavimo potrebne konstante za izračun F statistike in Fisherjeve porazdelitve  (q, n-k_u; alpha)
k_u <- 4
k_r <- 1
q <- k_u - k_r
n <- nrow(Z)
alpha <- 0.05

# poračunamo vrednost F statistike in Fisherjeve porazdelitve
F_stat <- ((VKR_r - VKR_u)/q)/(VKR_u/(n-k_u))
fisher <- qf(1-alpha, df1 = q, df2 = n - k_u)

# naredimo test
F_stat > fisher

# Zavrnemo ničelno hipotezo.

###############################################################################
# a) (ii) na 2. način - iz https://www.mattblackwell.org/files/teaching/ftests.pdf
n <- nrow(Z)
d <- ncol(Z)

L <- matrix(c(0,0,0,1,0,0,0,1,0,0,0,1), nrow=3, ncol=4)
c <- c(0,0,0)
H0 <- L %*% beta_hat_posploseni - c
q <- nrow(L)

varofdev <- solve(L %*% ztz_inv %*% t(L))
F_0 <- (t(H0) %*% varofdev %*% H0) / q
VKR <- t((X - Z %*% beta_hat_posploseni)) %*% (X - Z %*% beta_hat_posploseni)
# izračunami F statistiko
F_stat2 <- F_0/(VKR/(n-d))
# zelo podobna kot tista poračunana z algoritmom

#test - fisherjeva porazdelitev je enaka kot prej
F_stat2 > fisher
# Tudi po tem načinu hipotezo H0 zavrnemo.

###############################################################################

# a) (ii) na 3. način (matrično)

r <- 3

# nastavimo beto pod ničelno hipotezo
beta_hat_H <- c(mean(X), 0, 0, 0)

# poračunamo VKR z beta_hat_ (poljuben iz točke (i))
VKR <- t((X - Z %*% beta_hat_posploseni)) %*% (X - Z %*% beta_hat_posploseni)

# poračunamo še F statistike in Fisherjeve porazdelitve  (r, n-d; alpha)
F_stat3 <- ((t(Z %*% beta_hat_posploseni - Z %*% beta_hat_H) %*% (Z %*% beta_hat_posploseni - Z %*% beta_hat_H)) / r ) /
  (VKR/(n - d))
# F_stat3 je enaka kot F_stat2
fisher <- qf(1-alpha, df1 = r, df2 = n - d) 

#test - fisher še vedno isti
F_stat3 > fisher
# Tudi po tem načinu hipotezo H0 zavrnemo.

############################################################################

# (iii) Preizkusite domnevo Beta0=0 na standarden način

# R zelo prijazno sam od sebe poračuna t-statistiko v lm modelu.
# Treba je samo ustrezno primerjati:

coeff_beta_df <- data.frame(beta_hat_lm)
t_stat <- unlist(coeff_beta_df["(Intercept)", "t.value"])
student <- qt(1-alpha/2, df = n - d)
# testiramo hipotezo
abs(t_stat) > student

# Torej H0 zavrnemo

# a) (iii) 2. način - po pdfju
n <- nrow(Z)
d <- ncol(Z)

L <- matrix(c(1,0,0,0), nrow=1, ncol=4)
c <- 0
H0 <- L %*% beta_hat_posploseni - c
q <- nrow(L)

varofdev <- solve(L %*% ztz_inv %*% t(L))
F_0 <- (t(H0) %*% varofdev %*% H0) / q
VKR <- t((X - Z %*% beta_hat_posploseni)) %*% (X - Z %*% beta_hat_posploseni)

# Poračunamo F statistiko in fisherjevo porazdelitev (q, n-d, alpha) - fisher je tu drugačen
F_stat2 <- F_0/(VKR/(n-d))
fisher <- qf(1-alpha, df1 = q, df2 = n - d)

#test
F_stat2 > fisher
# Tudi po tem načinu hipoteze H0 ne zavrnemo.


# a)(iii) 3. način - matrično s t-testom
v <- c(1, 0, 0, 0)
t_stat2 <- abs(t(beta_hat_posploseni) %*% v / ( sqrt(t(ztz_inv %*% v) %*% v) %*% sqrt(VKR / (n-d))))
student <- qt(0.975, df = n - d)

t_stat2 > student
# torej hipoteze H0 ne zavrnemo

###############################################################################
# b) Predpostavimo linerani model brez prostega člena
# LDL = Beta_1 * TCH + Beta_2 * HDL + Beta_3 * TRI + epsilon

# (i) Ocenite Beta0, Beta1, Beta2, Beta3 po MNK - ponovimo postopek iz a), le brez prostega člena

# definiramo Z, X in epsilon = ei, z sigma = 1
sigma <- 1
ei <- rnorm(300, 0, sigma^2)

Z.2 <- as.matrix(podatki5[c(1:3)])

# prvi način - GS algoritem

#nastavimo mejo
eps <- 10^-10 
razcep <- grahm_schmidt(Z.2, eps)
S.2 <- razcep$S
P.2 <- razcep$P

norme_S.2 <- apply(S.2, 2, norma2)
J.2 <- diag(norme_S.2)

ztz.2 <- t(P.2) %*% J.2 %*% P.2
ztz_inv.2 <- solve(P.2) %*% J.2 %*% solve(t(P.2))

beta_hat_posploseni_2 <- ztz_inv.2 %*% t(Z.2) %*% X


# drugi način - navaden inverz

beta_hat_navadni_2 <- solve(t(Z.2) %*% Z.2) %*% t(Z.2) %*% X

# tretji način - z lm

lm.model.2 <- lm(data=podatki5, LDL ~ 0 + TCH + HDL + TRI) 
summary(lm.model.2)
beta_hat_lm_2 <- summary(lm.model.2)$coeff

###############################################################################

# (ii) Preizkusite domnevo (Beta1,Beta2,Beta3) = (1,-1,-0.45) na standarden način

# 1. način iz a) tu ne deluje
# po načinu iz pdf-ja

n <- nrow(Z.2)
d <- ncol(Z.2)

L.2 <- matrix(c(1,0,0,0,1,0,0,0,1), nrow=3, ncol=3)
c.2 <- c(1,-1,-0.45)
H02 <- L.2 %*% beta_hat_posploseni_2 - c.2
q <- nrow(L.2)

varofdev.2 <- solve(L.2 %*% ztz_inv.2 %*% t(L.2))
F_0 <- (t(H02) %*% varofdev.2 %*% H02) / q
VKR <- t((X - Z.2 %*% beta_hat_posploseni_2)) %*% (X - Z.2 %*% beta_hat_posploseni_2)

# poračunamo F statistiko in fisherjevo porazdelitev (q, n-d, alpha)
F_stat <- F_0/(VKR/(n-d))
fisher <- qf(1-alpha, df1 = q, df2 = n - d)

#test
F_stat > fisher
# Hipotezo H0 zavrnemo.


# b) (ii) na 2. način (matrično)
r <- 3

# nastavimo beto pod ničelno hipotezo
beta_hat_H <- c(1, -1, -0.45)

# VKR je enak kot v prejšnjem načinu
# poračunamo še vrednost  F statistike in Fisherjeve porazdelitve (r, n-d; alpha)
F_stat2 <- ((t(Z.2 %*% beta_hat_posploseni_2 - Z.2 %*% beta_hat_H) %*% 
               (Z.2 %*% beta_hat_posploseni_2 - Z.2 %*% beta_hat_H)) / r ) /
            (VKR/(n - d))
# F_stat2 je skoraj popolnoma enaka F_stat
fisher <- qf(1-alpha, df1 = r, df2 = n - d) 

#test
F_stat2 > fisher
# Tudi po tem načinu hipotezo H0 zavrnemo.

################################################################################

# b) (iii) Poišči območje zaupanja za (Beta1, Beta2, Beta3) stopnje zaupanja 0.95

# moramo narediti diagonalizacijo A = (Z^T)%*%Z.
# Na srečo ima R vgrajeno priročno funkcijo eigen prek katere lahko poiščemo
# lastne vektorje in lastne vrednosti matrike
A <- t(Z.2) %*% Z.2
# matrika lastnih vektorjev
Q <- eigen(A)$vectors
Q
# matrika lastnih vrednosti
lambda <- diag(eigen(A)$values)
lambda
# preverimo, če res deluje (mora biti enako A)
Q %*% lambda %*% solve(Q)
# res deluje, vse je OK

# izračunamo še ksi
ksi <- fisher * (d / (n - d)) * VKR
ksi

## POPRAVLJENO ######
# dolžine polosi
sqrt(ksi/lambda[1,1])
sqrt(ksi/lambda[2,2])
sqrt(ksi/lambda[3,3])

#############################################################################
# poskusimo narisati to območje zaupanja
library(plot3D)

# funkcija, ki spremeni koordinate v kartezične za elipsoid
spher2cart <- function(r,theta,phi) {
  
  x <- r * sin(theta) * cos(phi)
  y <- r * sin(theta) * sin(phi)
  z <- r * cos(theta)
  
  return(list(x = x, y = y, z = z))
}

# izračunamo dolžine osi ## POPRAVLJENO ##
a <- as.numeric(sqrt(ksi/lambda[1,1]))
b <- as.numeric(sqrt(ksi/lambda[2,2]))
c <- as.numeric(sqrt(ksi/lambda[3,3]))

# nastavimo parametre
theta <- seq(-pi/2, pi/2, length = 100)
phi <- seq(0, 2*pi, length = 100)
M <- mesh(theta, phi)
names(M) <- c("theta", "phi")

# parametrizacija za radij iz wikipedije https://en.wikipedia.org/wiki/Ellipsoid
R <- (a*b*c)/(sqrt(c^2*(b^2*cos(M$phi)^2 + a^2*sin(M$phi)^2)*cos(theta)^2 + a^2*b^2*sin(M$theta)^2))

cart <- spher2cart(R, M$theta, M$phi)

# narišemo 
par(mar = c(0, 0, 0, 0))
surf3D(cart$x, cart$y, cart$z, border = "black",
       colkey = FALSE, bty = "f",
       phi = 20, theta = 30)

# zgleda zanimivo :D

library(rgl)
plot3d(cart$x, cart$y, cart$z, size=3, 
       main = "Prikaz obmocja zaupanja stopnje zaupanja 0.95",
       xlab="TCH", ylab="HDL", zlab="TRI")

# funkcija za shranjevanje posnetka zaslona rgl
#rgl.snapshot('naloga5biii_obmocje_zaupanja2.png', fmt = 'png')
#################################################################################
# b) (iv)

t_bonferroni <- qt(1 - 0.05/(2*d), n - d)

# izračunamo delte za vsak Beta_j posebej
meja1 <- t_bonferroni * sqrt(ztz_inv.2[1, 1] * VKR /(n - d))
meja2 <- t_bonferroni * sqrt(ztz_inv.2[2, 2] * VKR /(n - d))
meja3 <- t_bonferroni * sqrt(ztz_inv.2[3, 3] * VKR /(n - d))

# zračunamo še nove intervale
IZ1 <- c(beta_hat_posploseni_2[1] - meja1, beta_hat_posploseni_2[1] + meja1)
IZ2 <- c(beta_hat_posploseni_2[2] - meja2, beta_hat_posploseni_2[2] + meja2)
IZ3 <- c(beta_hat_posploseni_2[3] - meja3, beta_hat_posploseni_2[3] + meja3)

# funkcija
IZ_bonferroni <- function(i, d, n, inverz, beta, VKR){
  t_bonferroni <- qt(1 - 0.05/(2*d), n - d)
  meja <- t_bonferroni * sqrt(inverz[i, i] * VKR /(n - d))
  IZ <- c(beta[i] - meja, beta[i] + meja)
  return(IZ)
}

IZ1 <- IZ_bonferroni(1,d,n,ztz_inv.2, beta_hat_posploseni_2, VKR)
IZ1
IZ2 <- IZ_bonferroni(2,d,n,ztz_inv.2, beta_hat_posploseni_2, VKR)
IZ2
IZ3 <- IZ_bonferroni(3,d,n,ztz_inv.2, beta_hat_posploseni_2, VKR)
IZ3

x_b <- seq(IZ1[1], IZ1[2], length = 50)
y_b <- seq(IZ2[1], IZ2[2], length = 50)
z_b <- seq(IZ3[1], IZ3[2], length = 50)

bonferronijevo_OZ <- expand.grid(x_b, y_b, z_b)
names(bonferronijevo_OZ) <- c("x_b", "y_b", "z_b")

plot3d(bonferronijevo_OZ$x_b,bonferronijevo_OZ$y_b,bonferronijevo_OZ$z_b, size=3,
       main = "Prikaz Bonferronijevega obmocja zaupanja stopnje zaupanja 0.95",
       xlab="TCH", ylab="HDL", zlab="TRI")

# Primerjava območij zaupanja
# primerjava volumna:
V_elipsoid = 4/3 * pi * a * b * c
V_elipsoid
V_kvader = (IZ1[2] - IZ1[1]) * (IZ2[2] - IZ2[1]) * (IZ3[2] - IZ3[1])
V_kvader
