###################################################################
# Predpostavimo model y = a1 * x + a2 * x^2 + a3 * x^3 + epsilon,
# kjer je epsilon porazdeljen N(0, sigma^2).
# 1. Oceni parametre a_i na dva načina:
#   a) preko formule (ZtZ)^-1 * Zt * X za primerni matriki X in Z
#   b) preko funkcije lm
# 2. Na graf nariši točke (x, y) ter krivuljo 
#   y = a1_hat * x + a2_hat * x^2 + a3_hat * x^3
# 3. Na podlagi ocen a_i_hat za parametre a_i napovej vrednosti y za x med
#    10 in 15. Napovedi dodaj na graf (z drugo barvo). 
###################################################################
set.seed(123)
n = 100
x = sort(runif(n, 0, 10))
epsilon = rnorm(100, 0, 1)
y = 0.3 * x - 0.5 * x^2 + 0.1 * x^3 + epsilon

#1.a)
X <- y
Z <- cbind(x, x^2, x^3)
aji <- solve(t(Z) %*% Z) %*% t(Z) %*% X

#b)
aji.lm <- lm(X~0 + Z)
summary(aji.lm)$coefficients
aji.hat <- aji.lm$coefficients


#2.)
plot(x,y)
y2 <- Z %*% aji.hat
lines(x,y2,col="red")

#3.)
x2 <- seq(10,15,0.01)
napoved <- aji.hat[1]*x2 + aji.hat[2]*x2^2 + aji.hat[3]*x2^3
plot(x,y, xlim=c(0,max(x2)), ylim=c(-1,max(napoved)))
lines(x,y2, col="red")
lines(x2,napoved, col="green")
legend("topleft", legend=c("podatki","ocena","napoved"), col=c("black", "red", "green"), lty=c(NA,1,1), pch=c(1,NA,NA))
title("Ocenjevanje in napovedovanje parametrov")
