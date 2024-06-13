#vaje 10
library(dplyr)

data <- read.csv("data.csv")
tabela <- table(data$izobrazba, data$soseska)
tabela

# v R:
chisq.test(tabela)
# če je p < 0.05 potem hipotezo ZAVRNEMO (to si želimo)
qchisq(0.95, 6)
#zavrnemo, ker po chisq.test X^2 > qchisq


uji <- rowSums(tabela)
vji <- colSums(tabela)
n <- sum(vji)


vrednost = function(i, j){
  T_hat <- uji[i]*vji[j]/n
  (tabela[i, j] - T_hat)^2/T_hat
}

vrednosti_par <- expand.grid(1:nrow(tabela), 1:ncol(tabela))

Map(function(i,j) vrednost(i,j), vrednosti_par[1] %>% unlist(), vrednosti_par[2] %>% unlist()) %>% 
  unlist() %>% sum()


# T2 - delimo z opaženo frekvenco namesto s pričakovano
vrednost2 = function(i, j){
  T_hat <- uji[i]*vji[j]/n
  (tabela[i, j] - T_hat)^2/tabela[i, j]
}
Map(function(i,j) vrednost2(i,j), vrednosti_par[1] %>% unlist(), vrednosti_par[2] %>% unlist()) %>% 
  unlist() %>% sum()


# razmerje verjetij
vrednost3 = function(i, j){
  T_hat <- uji[i]*vji[j]/n
  tabela[i,j] * log(T_hat / tabela[i, j])
}
tL = Map(function(i,j) vrednost3(i,j), vrednosti_par[1] %>% unlist(), vrednosti_par[2] %>% unlist()) %>% 
  unlist() %>% sum()
-2*tL

###############################################################

#preizkus homogenosti - vaje 11
druzine <- c(5, 15, 35, 17, 28)
samski <- c(45, 65, 37, 46, 7)

rezultati <- as.data.frame(rbind(druzine, samski))
colnames(rezultati) <- c("sportni", "limuzina", "kombi limuzina", "tovornjak", "SUV")
chisq.test(rezultati)
# zavrnemo hipotezo homogenosti -> porazdelitev avtov ni enaka
qchisq(0.95, 4)

# poračunajmo sami
tabela = rezultati

uji <- rowSums(tabela)
vji <- colSums(tabela)
n <- sum(vji)


vrednost = function(i, j){
  T_hat <- uji[i]*vji[j]/n
  (tabela[i, j] - T_hat)^2/T_hat
}

vrednosti_par <- expand.grid(1:nrow(tabela), 1:ncol(tabela))

Map(function(i,j) vrednost(i,j), vrednosti_par[1] %>% unlist(), vrednosti_par[2] %>% unlist()) %>% 
  unlist() %>% sum()
