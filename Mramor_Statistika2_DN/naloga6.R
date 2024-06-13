# 6. naloga

podatki6 <- read.delim("podatki_6.txt", header=FALSE)

colnames(podatki6) <- c("weight", "mpg", "foreign")

# b) Ocenite parametre Beta0, Beta1 in Beta2 po MNV

beta <- rep(0,3)

Z <- podatki6[, -3]
Z <- as.matrix(cbind(1, Z))
X <- podatki6[,3]

k <- 1
newton_raphson <- function(beta, Z, X){
  while(k > 0.000001){
    z_beta <- apply(Z, 1, function(y) y %*% beta)
    p <- 1-exp(-exp(z_beta))
    A <- (X/p - 1) * exp(z_beta)
    
    gradient <- t(Z) %*% A

    H <- t(Z) %*% diag(exp(z_beta)*(X/p - X/p^2 * (1-p) * exp(z_beta) - 1)) %*% Z
    
    
    k <- solve(H) %*% gradient
    beta <- beta - k
    k <- max(abs(k))
  }
  return(beta)
}

bete_MNV <- newton_raphson(beta, Z, X)
bete_MNV

# vemo da je naša funkcija verjetnosti complementary log-log 
# zato lahko ocenimo bete z glm za primerjavo
model <- glm(formula = foreign ~ weight + mpg, data = podatki6, family = binomial(link=cloglog))
beta_prave <- model$coefficients
beta_prave
# izgleda približno podobno - verjetno je naš približek primeren.

################################################################################

# c) Fisherjeva informacijska matrika

z_beta <- apply(Z, 1, function(y) y%*% beta_prave)
p <- 1-exp(-exp(z_beta))
E_H <- t(Z) %*% diag(- exp(z_beta)^2 * (1-p)/p) %*% Z
n <- length(X)

FI = - E_H/n
FI

################################################################################

# d) Izračunajte standardne napake za parametre Beta_0, Beta_1, Beta_2

#stanardni odklon
var_beta.hat <- solve(FI)/n
sqrt(var_beta.hat[1,1])
sqrt(var_beta.hat[2,2])
sqrt(var_beta.hat[3,3])

summary(model)$coefficients

# enake vrednosti kot pri glm modelu, torej je vse v redu.

###############################################################################

# e) Na stadndarden način preizkusite H0: Beta_1=Beta_2=0

# testiranje hipotez 
z_beta <- apply(Z, 1, function(y) y%*% bete_MNV)
l_full = sum(X*log(1-exp(-exp(z_beta))) - (1-X)*exp(z_beta))

# pod H0
beta0 <- log(-log(1-mean(X)))
glm(formula = foreign ~ 1, data = podatki6, family = binomial(link=cloglog))

beta2 <- c(beta0, 0, 0)
z_beta2 <- apply(Z, 1, function(y) y%*% beta2)
l_reduced = sum(X*log(1-exp(-exp(z_beta2))) - (1-X)*exp(z_beta2))

# izračunamo testno statistiko po metodi razmerja verjetij
lambda <- -2*(l_reduced - l_full)
lambda
hi2 <- qchisq(0.95, 2)
hi2

lambda > hi2
# Ker je testna statistika > chisq, hipotezo 0 zavrnemo.

model0 <- glm(formula = foreign ~ weight + mpg, data = podatki6, family = binomial(link=cloglog))
model1 <- glm(formula = foreign ~ 1, data = podatki6, family = binomial(link=cloglog))

anova(model1, model0)

# deviance med modeloma se ujema z našo testno statistiko - vse je OK

