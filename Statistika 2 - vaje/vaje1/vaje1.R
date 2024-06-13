# 1. povp volumen kave 

mu = 40
sd = 2
n = 10000

set.seed(123)
data1 <- rnorm(n, mu, sd)
hist(data1)

#MM - metoda momentov
mu1 = mean(data1)
sd1 = mean(data1^2) - mean(data1)^2

#MNV - metoda največjega verjetja
mu2 = mean(data1)
sd2 = mean((data1-mu2)^2)

var(data1)


#2. eksponentne družine 
a=2
lambda=3
data2 = rgamma(n,a,lambda)
