#vaje12

before = c(25,28,32,18,15,18,25,31)
after = c(37,26,47,25,29,22,19,32)

data = data.frame(before,after)
#razlike: prej-potem
data$diff = data$before - data$after
#absolutna vrednost razlike
data$abs = abs(data$diff)
# rank
data$rank = rank(data$abs)

#vsota pozitivnih rankov
library(tidyverse)
data %>% filter(diff > 0) %>% summarise(s = sum(rank))


# p-value iz predavanj
p = 14/256
p
# p-value v R-ju
w = wilcox.test(before, after, paired=TRUE,
            alternative = "less")
w$p.value
#p == w$p.value


# recursive formula

u <- function(n, t){
    if(n == 0 & t == 0){
      1
    }
    else if(n == 0 & t != 0){
      0
    }
    else if(t < 0 | t > n*(n+1)/2){
      0
    }
    else{
      u(n-1, t) + u(n-1, t-n)
    }
}

# P(RS <= 6 w mu = 0)
n = nrow(data)
sapply(0:6, function(x) u(n,x)) %>% sum() /2^n

vsota <- 0
for (i in (0:6)){
  vsota <- vsota + u(8,i)
}
vsota


# critical value
# C max, s.t. P(RS <= C) <= alpha
sapply(0:5, function(x) u(n,x)) %>% sum() /2^n
# C = 5


alpha = 0.05
for (i in 0:(n*(n+1)/2)){
  p = sapply(0:i, function(x) u(n,x)) %>% sum() /2^n
  if (p >= alpha){
    C = i-1
    break #našli smo C
  }
}
C


##########################################################
data = read.csv("blood_pressure.csv")

#razlike: prej-potem
data$diff = data$pritisk_pred - data$pritisk_po
#absolutna vrednost razlike
data$abs = abs(data$diff)
# rank
data$rank = rank(data$abs)

#vsota pozitivnih rankov
data %>% filter(diff > 0) %>% summarise(s = sum(rank))

# we test mu <= 0 (radi bi zavrnili to hipotezo)
wilcox.test(data$pritisk_pred, data$pritisk_po,
            paired=TRUE, alternative = "greater")
# zavrnemo hipotezo mu <= 0, ker je p-value < 0.05

# test normality of data
# H0: data is normal
shapiro.test(data$pritisk_pred)
shapiro.test(data$pritisk_po)

# ne moremo zavrnit hipoteze, da so podatki normalno porazdeljeni
# uporabimo t-test
# we test mu <= 0 (radi bi zavrnili to hipotezo)
t.test(data$pritisk_pred, data$pritisk_po,
            paired=TRUE, alternative = "greater")
# tudi t nam zavrne hipotezo, da je mu <= 0

# sign test
n = nrow(data)
diff_plus = data %>% filter(diff>0) %>% nrow()
# testiramo p+ <= p- (želimo zavrnit) <=> q < 1/2
binom.test(x=diff_plus, n=n, p=0.5,
           alternative = "greater")

#p-value
#P(Bin(n, 0.5) >= diff_plus)
sapply(diff_plus:n, function(x) dbinom(x, n, 1/2)) %>% sum()


# asymptotic test (Z-test)
RS = data %>% filter(diff > 0) %>% summarise(s = sum(rank))
#
Z = (RS - n*(n+1)/4)/sqrt(n*(n+1)*(2*n+1)/24)
Z = as.numeric(Z)

wilcox.test(data$pritisk_pred, data$pritisk_po,
            paired=TRUE, alternative = "greater")
qnorm(0.95)

Z > qnorm(0.95)
# zavrnemo ničelno hipotezo

# p value
pnorm(-Z)
1-pnorm(Z)
wilcox.test(data$pritisk_pred, data$pritisk_po,
            paired=TRUE, alternative = "greater",
            exact = FALSE, correct = FALSE)
##########################################
# 13. vaje
data = data.frame(
  participant = 1:12,
  product = c(rep("X", 6), rep("Y", 6)),
  value = c(3,4,2,6,2.5,5,9,7,5.5,10,6.5,8)
)
data$rank = rank(data$value)
library(tidyverse)
RS = data %>% filter(product == "X") %>% summarise(s = sum(rank))

# asymptotic test
m= 6
n= 6
exp_value = m * (m+n+1)/2
exp_value
var = m * n * (m+n+1)/12
var

# Z test
Z = (RS - exp_value)/sqrt(var)
Z
alpha = 0.05
z_alpha_2 = qnorm(1- alpha/2)
# zavrnemo, če je |Z| > z_alpha_2
abs(Z) > z_alpha_2

# zavrnemo H = porazdelitev vrednosti X in Y je enaka

# p value?
# verjetnost, da ima Z vrednost ekstremnejšejšo od -2.722179
pnorm(as.numeric(Z)) * 2


x = data %>% filter(product  == "X") %>% select(value) %>% unlist()
y = data %>% filter(product  == "Y") %>% select(value) %>% unlist()

#exact test
wilcox.test(x, y, alternative = "two.sided")

#asymptotic test
wilcox.test(x, y, alternative = "two.sided",
            exact = FALSE, correct = FALSE)

###############################################
# enak example but with ties

##########################################
# 13. vaje
data = data.frame(
  participant = 1:12,
  product = c(rep("X", 6), rep("Y", 6)),
  value = c(3,4,2,6,2,5,9,7,5,10,6,8)
)
data$rank = rank(data$value)
library(tidyverse)
RS = data %>% filter(product == "X") %>% summarise(s = sum(rank))

# asymptotic test
m= 6
n= 6

# correction for ties
exp_value = m * (m+n+1)/2
exp_value

rank_number = table(data$rank)
rank_ties = rank_number[rank_number > 1]
vsota = sum(rank_ties^3 - rank_ties)

correction = m * n * sum(rank_number^3 - rank_number)/(12*(m+n)*(m+n-1))

var = m * n * (m+n+1)/12 - correction
var

# Z test
Z = (RS - exp_value)/sqrt(var)
Z
alpha = 0.05
z_alpha_2 = qnorm(1- alpha/2)
# zavrnemo, če je |Z| > z_alpha_2
abs(Z) > z_alpha_2

# zavrnemo H = porazdelitev vrednosti X in Y je enaka

# p value?
# verjetnost, da ima Z vrednost ekstremnejšejšo od -2.722179
pnorm(as.numeric(Z)) * 2


x = data %>% filter(product  == "X") %>% select(value) %>% unlist()
y = data %>% filter(product  == "Y") %>% select(value) %>% unlist()

#exact test - warning
wilcox.test(x, y, alternative = "two.sided")

#asymptotic test
wilcox.test(x, y, alternative = "two.sided",
            exact = FALSE)

# same
wilcox.test(x, y, alternative = "two.sided",
            exact = FALSE, correct = FALSE)
pnorm(as.numeric(Z)) * 2
