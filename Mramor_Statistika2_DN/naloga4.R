# 4. naloga

library(dplyr)
library(tidyverse)
library(partitions)


p <- c(1/10, 1/10, 2/5, 1/5, 1/5)
m <- 5
n <- c(30, 50, 70, 90)

# najprej potrebujemo vse možne izide multinomske porazdelitve za vsak n pri danem m

# 1. poskus - časovno zelo zahtevno - ne deluje za n = 70 in n = 90 :(
kombinacije <- function(n){
  tabela <- as.matrix(expand.grid(0:n, 0:n, 0:n, 0:n, 0:n))
  tabela <- tabela[rowSums(tabela) == n, ]
  return(tabela)
}

#comb30 <- kombinacije(30)
#comb50 <- kombinacije(50)
#comb70 <- kombinacije(70)
#comb90 <- kombinacije(90)

# 2. poskus - s compositions

kombinacije2 <- function(n,m){
  tabela <- t(as.matrix(compositions(n,m)))
  return(tabela)
}


comb30 <- as.data.frame(kombinacije2(n[1],m))
comb50 <- as.data.frame(kombinacije2(n[2],m))
comb70 <- as.data.frame(kombinacije2(n[3],m))
comb90 <- as.data.frame(kombinacije2(n[4],m))

# to dela odlično!

##############################################################################

# a) velikost preizkusa domneve na podlagi razmerja verjetij

# testna verzija za n = 30
chi <- qchisq(0.95, m-1)

comb30 <- comb30 %>% rowwise() %>% mutate(verjetnost = dmultinom(c(V1,V2,V3,V4,V5), prob = p)) %>% 
          mutate(lambda = prod( (n[1]*p/c(V1,V2,V3,V4,V5))^c(V1,V2,V3,V4,V5))) %>% 
          mutate(statistika = -2*log(lambda)) %>% mutate(test = statistika > chi)

zavrnemo <- comb30[comb30$test == TRUE,]
sum(zavrnemo$verjetnost)

##################################################################################

# sedaj sestavimo funkcijo za računanje v splošnem

velikost_preizkusa <- function(n, p, m){
  tabela <- as.data.frame(kombinacije2(n,m))
  chi <- qchisq(0.95, m-1)
  
  tabela <- tabela %>% rowwise() %>% 
    mutate(verjetnost = dmultinom(c(V1,V2,V3,V4,V5), prob = p)) %>% 
    mutate(lambda = prod( (n*p/c(V1,V2,V3,V4,V5))^c(V1,V2,V3,V4,V5))) %>% 
    mutate(statistika = -2*log(lambda)) %>% mutate(test = statistika > chi)
  
  zavrnemo <- tabela[tabela$test == TRUE,]
  return(sum(zavrnemo$verjetnost))
}

# poračunajmo velikosti preizkusa za vse n
# pri velikih n je treba biti potrpežljiv, saj je problem računsko precej zahteven
# in lahko traja nekoliko dlje, da program najde velikost preizkusa.
velikost_preizkusa30 <- velikost_preizkusa(30, p, m)
velikost_preizkusa30
velikost_preizkusa50 <- velikost_preizkusa(50, p, m)
velikost_preizkusa50
velikost_preizkusa70 <- velikost_preizkusa(70, p, m)
velikost_preizkusa70
velikost_preizkusa90 <- velikost_preizkusa(90, p, m)
velikost_preizkusa90

#############################################################################

# b) velikost preizkusa domneve na podlagi Pearsonove statistike

#postopamo enako kot prej - v funkciji spremenimo le testno statistiko

velikost_preizkusa_Pearson <- function(n, p, m){
  tabela <- as.data.frame(kombinacije2(n,m))
  chi <- qchisq(0.95, m-1)
  
  tabela <- tabela %>% rowwise() %>% 
    mutate(verjetnost = dmultinom(c(V1,V2,V3,V4,V5), prob = p)) %>% 
    mutate(statistika = sum( (c(V1,V2,V3,V4,V5) - n*p)^2 / (n*p) )) %>% 
    mutate(test = statistika > chi)
  
  zavrnemo <- tabela[tabela$test == TRUE,]
  return(sum(zavrnemo$verjetnost))
}


# poračunajmo velikosti preizkusa za vse n
velikost_preizkusa_Pearson30 <- velikost_preizkusa_Pearson(30, p, m)
velikost_preizkusa_Pearson30
velikost_preizkusa_Pearson50 <- velikost_preizkusa_Pearson(50, p, m)
velikost_preizkusa_Pearson50
velikost_preizkusa_Pearson70 <- velikost_preizkusa_Pearson(70, p, m)
velikost_preizkusa_Pearson70
velikost_preizkusa_Pearson90 <- velikost_preizkusa_Pearson(90, p, m)
velikost_preizkusa_Pearson90

# Rezlutati delujejo primerno