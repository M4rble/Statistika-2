# vaje 11

# TEST z ZNAKI
n = 10
t0 = 0
binom.test(8, n-t0, p = 0.5, alternative = "two.sided")
# ne zavrnemo hipoteze, da je p+ = p-

# iščemo C
alpha = 0.05
C = 7
# ta verjetnost se veča z večjim C
pbinom(C, n-t0, 1/2)
# iskanje na roke? neeee


meja = 1- alpha/2
for (i in 0:(n-t0)){
  p = pbinom(i, n-t0, 1/2)
  if (p > meja){
    C = i
    break
  }
}
C

Tplus = 8 # števlo pozitivnih razlik v dolžini nog

# zavrnemo hipotezo p+ = p-, če
Tplus < n -t0 - C
# ali 
Tplus > C


#C = 0
#while (C < 11){
#  p = pbinom(C, n-t0, 1/2)
#  if (p > meja){
#    break
#  }
#  C = C + 1
#}
#C

