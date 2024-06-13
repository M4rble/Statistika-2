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
n = 100
x = runif(n, 0, 10)
epsilon = rnorm(100, 0, 8)
y = 0.3 * x + 1.5 * x^2 + 0.1 * x^3 + epsilon

