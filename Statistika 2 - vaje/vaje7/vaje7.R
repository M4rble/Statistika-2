#vaje 7

podatki <- read.csv("logit.csv")
n <- nrow(podatki)

#ocenjevanje z glm
m1 <- glm(formula = x ~ z1 + z2, data = podatki, family = binomial(link=logit))
m1$coefficients

#ocenjevanje z newton-raphsonom
beta_0 = rep(0,3)

Z <- podatki[, -1]
Z <- as.matrix(cbind(1, Z))
X <- podatki[,1]

beta_minus1 = rep(1,3)

newton_raphson <- function(beta_0, Z, X){
  beta_t1 <- rep(1,3)
  beta_t <- beta_0
  while(abs(beta_t - beta_t1) > rep(10^-5,length(beta_0))){
    p <- rep(0, 1000)
    for (i in (1:1000)){
      p[i] <- exp(Z[i,] %*% beta_t)/(1+exp(Z[i,] %*% beta_t))
    }
    

    grad_beta <- t(Z) %*% (X - p)
    
    p_diag <- matrix(0, nrow=1000, ncol=1000)
    for(i in (1:1000)){
      p_diag[i,i] <- p[i]*(1-p[i])
    }
    H_beta <- - t(Z) %*% p_diag %*% Z
    
    beta_t1 <- beta_t
    beta_t <- beta_t - solve(H_beta) %*% grad_beta
  }
  return(beta_t)
}

#newton_raphson(beta_0, Z, X)



logit_inv <- function(x){
  exp(x)/(1+ exp(x))
}
  

k <- 1
beta <- rep(0,3)
while(k > 0.00001){
  z_beta <- apply(Z, 1, function(y) y%*% beta)
  p_ji <- logit_inv(z_beta)
  H <- - t(Z) %*% diag(p_ji * (1-p_ji)) %*% Z
  grad <- t(Z) %*% (X - p_ji)
  
  k <- solve(H) %*% grad
  beta <- beta - k
  k <- max(abs(k))
}
beta

# fisherjeva informacijska matrika
z_beta <- apply(Z, 1, function(y) y%*% beta)
p_ji <- logit_inv(z_beta)
H <- - t(Z) %*% diag(p_ji * (1-p_ji)) %*% Z

FI = - H/n

#stanardni odklon
var_beta.hat <- solve(FI)/n
sqrt(var_beta.hat[1,1])
sqrt(var_beta.hat[2,2])
sqrt(var_beta.hat[3,3])

summary(m1)$coefficients

# testiranje hipotez 
z_beta <- apply(Z, 1, function(y) y%*% beta)
l_full = sum(-log(1 + exp(z_beta)) + X*z_beta)

# pod H0
beta0 <- log(mean(X)/(1-mean(X)))
glm(formula = x ~ 1, data = podatki, family = binomial(link=logit))

beta2 <- c(beta0, 0, 0)
z_beta2 <- apply(Z, 1, function(y) y%*% beta2)
l_reduced = sum(-log(1 + exp(z_beta2)) + X*z_beta2)

# nekaj ni ok, to bi moglo prit enako kot deviance po anovi
-2*(l_reduced - l_full)
# je vse OK?

model0 <- glm(formula = x ~ z1 + z2, data = podatki, family = binomial(link=logit))
model1 <- glm(formula = x ~ 1, data = podatki, family = binomial(link=logit))

anova(model1, model0)



######################################

# vaje 8
#ocenjevanje z glm
m2 <- glm(formula = x ~ z1 + z2, data = podatki, family = binomial(link=probit))
m2$coefficients

beta = rep(0,3)
k <- 1
q <- 2*X - 1
while(k > 0.00001){
  z_beta <- apply(Z, 1, function(y) y %*% beta)
  q.z_beta <- q*z_beta
  
  gostota <- dnorm(q.z_beta)
  porazdelitev = pnorm(q.z_beta)
  lambda <- (gostota/porazdelitev) * q
  
  W <- diag(lambda * (z_beta + lambda))
  gradient <- t(Z) %*% lambda
  H <- - t(Z) %*% W %*% Z
  
  k <- solve(H) %*% gradient
  beta <- beta - k
  k <- max(abs(k))
}
beta


#vaje9 - FI in std odklon za probit

# fisherjeva informacijska matrika
z_beta <- apply(Z, 1, function(y) y%*% beta) 

W <- diag(dnorm(z_beta)^2 / (pnorm(z_beta)*(1-pnorm(z_beta))))
# - E[H(logL)]
E_H_minus <- - t(Z) %*% W %*% Z

FI = - E_H_minus/n

#stanardni odklon
var_beta.hat <- solve(FI)/n
sqrt(var_beta.hat[1,1])
sqrt(var_beta.hat[2,2])
sqrt(var_beta.hat[3,3])

summary(m2)$coefficients

#vaje 10 - testiranje hipotez

# full model
z_beta <- apply(Z, 1, function(y) y%*% beta)
l_full = sum(log(pnorm(q*z_beta)))

# CNV za beta0 pod H0
beta0 = 0

k = 1
while(k > 0.00001){
  lambda = dnorm(q * beta0)/pnorm(q * beta0) * q
  prvi_odvod = sum(lambda)
  drugi_odvod = - sum(lambda * (beta0 + lambda))
  k =  prvi_odvod/drugi_odvod
  beta0 = beta0 - k
  k = abs(k)
}
beta0

# beta pod H0
beta2 = c(beta0, 0, 0)

# pod H0 - reduced model
z_beta2 <- apply(Z, 1, function(y) y%*% beta2)
l_reduced = sum(log(pnorm(q*z_beta2)))

#vrednost testne statistike
-2 * (l_reduced - l_full)

qchisq(0.95, 2)

# H0 zavrnemo, ker je 503.3526 > od 5.991465 (-2ln(lambda) > hi^2_(2;0.95))

model0 = glm(formula = x ~ z1 + z2, data = podatki, family = binomial(link=probit))

model1 = glm(formula = x ~ 1, data = podatki, family = binomial(link=probit))
model1$coefficients
beta0


anova(model1, model0)
