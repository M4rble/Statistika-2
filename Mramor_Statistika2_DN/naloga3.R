# 3. naloga

library(ggplot2)
library(dplyr)


# b) Izračunaj gamma_1 in gamma_2
n <- 23
p0 <- 6/10
C1 <- 9
C2 <- 18
alpha <- 0.05


# funkcije za odvod in vsoto odvodov
odvod <- function(i, n, p){
  choose(n, i) * p^(i - 1) * (1 - p)^(n - i-1) * (i - n*p)
}

vsota_odv <- function(m, M, n, p){
  vsota <- 0
  for(i in m:M){
    vsota <- vsota + odvod(i, n, p)
  }
  return(vsota)
}

# reševali bomo matrično
# poračunamo elemente vsake matrike posebej

# elementi matrik A in B
# pbinom = P(X <= k) 
# dbinom = P(X = k)
a1 <- pbinom(C1, n, p0) - dbinom(C1, n, p0)
a2 <- 1 - pbinom(C2, n, p0)
a3 <- dbinom(C1, n, p0)
a4 <- dbinom(C2, n, p0)

b1 <- vsota_odv(0, C1 - 1, n, p0)
b2 <- vsota_odv(C2 + 1, n, n, p0)
b3 <- odvod(C1, n, p0)
b4 <- odvod(C2, n, p0)

# matrika A
A <- matrix(c(a3, b3, a4, b4), ncol = 2, nrow = 2)

# Matrika B
B <- c(alpha - a1 - a2, - b1 - b2)

gamma <- solve(A, B)
gamma

##############################################################################

# c) Narišite graf pokritosti in ocenite koeficient zaupanja

# interval zaupanja iz a)
S <- c(0.001, 0.002, 0.004, 0.016, 0.038, 0.066, 0.096, 0.128, 0.158, 0.192, 0.228,
       0.266, 0.304, 0.344, 0.384, 0.426, 0.470, 0.516, 0.564, 0.614, 0.668, 0.724, 0.786, 0.858)
Z <- c(0.144, 0.216, 0.278, 0.334, 0.386, 0.436, 0.484, 0.530, 0.574, 0.616, 0.656,
       0.698, 0.736, 0.772, 0.808, 0.842, 0.872, 0.906, 0.936, 0.962, 0.984, 0.996, 0.998, 0.999)

# razdelimo interval verjetnosti
N <- 9990

pokritost <- numeric(N)
verjetnosti <- seq(0.001,0.999,length=N)

for(i in 1:N){
  for(j in 1:length(S)){
    if(S[j] <= verjetnosti[i] & verjetnosti[i] <= Z[j]){
      pokritost[i] <- pokritost[i] + dbinom(j-1, n, verjetnosti[i])
      #print(pokritost)
    }}}

koef_zaupanja <- min(pokritost)
koef_zaupanja

graf <- ggplot(data.frame(verjetnosti, pokritost), 
               aes(x = verjetnosti, y = pokritost)) + 
  geom_point(size=0.5, shape=16, color = 'darkblue') +
  xlab("Verjetnosti") + ylab("Pokritost") + ggtitle('Graf pokritosti') + 
  theme_bw()
graf


##############################################################################

# d) za n1 = n+10 poiščite C1, C2, gamma1, gamma2

n1 <- n + 10

# definirajmo funkcijo, ki bo poiskala iskane konstante C1, C2, gamma1, gamma1
poisci_konstante <- function(n1, p0, alpha){
  
  # Gremo čez vse možne vrednosti za C1 in C2 ter gledamo možne kombinacije
  for(C1 in 0:n1){
    for(C2 in (C1+1):n1){
      
      # A konstruiramo enako kot v a)
      a3 <- dbinom(C1, n1, p0)
      a4 <- dbinom(C2, n1, p0)
      b3 <- odvod(C1, n1, p0)
      b4 <- odvod(C2, n1, p0)
      
      A <- matrix(c(a3, b3, a4, b4), ncol = 2, nrow = 2)
      
      # B konstruiramo podobno kot v B, pazimo na robne primere
      # robni primeri vsaj eden izmed C1,C2 je enak bodisi 0 ali n1
      if(C1 == 0 & C2 == n1){
        a1 <- 0
        a2 <- 0
        b1 <- 0
        b2 <- 0}
      
      if(C1 == 0 & C2 != n1){
        a1 <- 0
        a2 <- 1 - pbinom(C2, n1, p0)
        b1 <- 0
        b2 <- vsota_odv(C2 + 1, n1, n1, p0)}
      
      if(C1 != 0 & C2 == n1){
        a1 <- pbinom(C1, n1, p0) - dbinom(C1, n1, p0)
        a2 <- 0
        b1 <- vsota_odv(0, C1 - 1, n1, p0)
        b2 <- 0}
      
      if(C1 != 0 & C2 != n1){
        # enako kot v b)
        a1 <- pbinom(C1, n1, p0) - dbinom(C1, n1, p0)
        a2 <- 1 - pbinom(C2, n1, p0)
        b1 <- vsota_odv(0, C1 - 1, n1, p0)
        b2 <- vsota_odv(C2 + 1, n1, n1, p0)}
      
      B <- c(alpha - a1 - a2, - b1 - b2)
      
      tryCatch(gamma <- solve(A,B))
      
      # preverimo pogoj za gamma
      if(gamma[1] >= 0 & gamma[2] >= 0 & gamma[1] <= 1 & gamma[2] <= 1){
        return(c(C1, C2, gamma[1], gamma[2]))
      }
    }
  }
}

poisci_konstante(n1, 6/10, 0.05)


###############################################################################

# e) 

M <- 999
p_vec <- seq(0.001,0.999,length=M)

# definiramo funkciji za iskanje C1 in C2

poisci_C1 <- function(N, p, alpha){
  for (C in 0 : N) {
    if(pbinom(C, N, p) >= alpha){
      return(C)
    }
  }}

poisci_C2 <- function(N, p, alpha){
  for (C in 0 : N) {
    if(1-pbinom(C, N, p) <= alpha){
      return(C)
    }
  }}

C_sp <- lapply(p_vec, poisci_C1, N=n1, alpha = 0.05/2) %>% unlist()
table(C_sp)

C_zg <- lapply(p_vec, poisci_C2, N=n1, alpha = 0.05/2) %>% unlist()
table(C_zg)

col <- rep(1:3, M/3)

df_C <- as.data.frame(cbind(C_sp, C_zg, p_vec, col))
df_C$col <- factor(df_C$col)

diagram <- ggplot(df_C) + geom_point(aes(x=p_vec, y=C_sp, colour = col)) + 
      geom_point(aes(x=p_vec, y=C_zg, colour = col)) + 
      scale_color_manual(values = c("darkblue", "darkgreen", "darkred")) + 
      theme_bw() + theme(legend.position = "none") + 
      xlab("Vrednosti p") + ylab("Vrednosti C") + ggtitle("Diagram kot v a) za n=33") +
      scale_x_continuous(breaks = seq(0,1,0.1)) + 
      scale_y_continuous(breaks = seq(0,33,1))
diagram
