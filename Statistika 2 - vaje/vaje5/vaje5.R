# Stastistika - Vaje 5
# Ocenjenvane funkcije gamma

#simulacija
a <- 7
b <- 0.5

vzorec <- rgamma(10000,a,b)

# Metoda momentov 

m1 <- mean(vzorec)
m2 <- mean(vzorec^2)

b.mm <- m1/(m2-m1^2)
a.mm <- m1^2/(m2-m1^2)


# Metotda največjega verjetja

g <- function(x) digamma(x) - log(x) - mean(log(vzorec)) + log(mean(vzorec))

g_odvod <- function(x) trigamma(x) - 1/x

a0 <- a.mm
for(i in 1:1000){
  a0 <- a0 - g(a0)/g_odvod(a0)
}

g(a0)

#uniroot
uniroot(g, c(0.00001, a.mm + 4))$root

b.mle <- a0/mean(vzorec)
b.mle
