# 2. naloga

# potrebne knjižnice
library(ggplot2)
library(dplyr)
library(GenBinomApps)

# uvozimo podatke
podatki2 <- read.csv("podatki_2.txt", header=FALSE)

# a) kompletna zadostna statistika

T_2 <- sum(podatki2$V1)

# c)
# nastavimo vrednosti
n <- 100
m <- nrow(podatki2)
N <- n*m
theta0 <- 0.2
alpha = 0.05

# sestavimo funkcijo - vzamemo C za katerega bo verj. porazdelitev <= alpha = 0.05
poisci_C <- function(N, theta0, alpha){
  for (C in 0 : 30000) {
    if(1-pbinom(C, N, theta0) <= alpha){
      return(C)
    }
  }}

C = poisci_C(N, theta0, alpha)
C

# izračunamo še gammo  
alpha_izr <- 1 - pbinom(C, N, theta0)
gamma <- (alpha - alpha_izr)/dbinom(C, N, theta0)
gamma

# preverimpo, da je st. znač res 0.05
#1 - pbinom(C, N, theta0) + gamma*dbinom(C, N, theta0)

########################################################################
# primerjamo z izračunom C po aproksimaciji z normalno porazdelitvijo

inv_norm <- qnorm(1-alpha)
C_2 <- inv_norm * sqrt(N * theta0 * (1-theta0)) + N*theta0
C_2 <- round(C_2)
C_2

gamma_2 <- (alpha - 1 + pbinom(C_2, N, theta0))/dbinom(C_2, N, theta0)
gamma_2
# izračunani gamma seveda pride enak kot po prvem postopku.
##########################################################################

# narišimo funkcijo moči

theta <- seq(0, 1, by = 0.0001)
vrednosti <- lapply(theta, function(theta) 1-pbinom(C, N , theta) + gamma * dbinom(C, N , theta)) %>% unlist()

graf_moci <- ggplot(data.frame(theta, vrednosti), aes(theta, vrednosti)) + geom_line( color = 'darkblue') +
  xlab("Theta") + ylab("Moč") + ggtitle('Graf funkcije moči') + 
  geom_vline(aes(xintercept = theta0), color ="red", lty = "dashed") +
  geom_hline(aes(yintercept = alpha), color = "black", lty = "dashed") + theme_bw() +
  annotate(geom="text", x=0.02, y=0.08, label="alpha = 0.05", col="black") + 
  annotate(geom="text", x=0.13, y=1, label="theta0 = 0.2", col="red") 

graf_moci

##########################################################################

# d) Interval zaupanja stopnje zaupanja 0.95

poisci_theta <- function(T, N, alpha){
  for (theta in seq(0,1,0.000001)) {
    if(pbinom(T, N, theta) <= 1-alpha){
      theta_sp = theta
      break()
    }}
  for (theta in seq(1,0,-0.000001)) {
    if(pbinom(T, N, theta) >= alpha){
        theta_zg = theta
        break()
    }}
  return(c(theta_sp, theta_zg))
}

poisci_theta(T_2, N, 0.475)


# Z vgrajeno metodo clopper-pearsonov IZ:

clopper.pearson.ci(T_2, N, alpha = 0.95, CI = "two.sided")

# Naš izračun se precej dobro ujema z vgrajeno Clopper-Pearsonovo metodo.
