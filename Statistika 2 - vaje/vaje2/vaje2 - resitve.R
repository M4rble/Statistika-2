###################################################################
# Splošni linearni model je oblike X = Z * beta + epsilon
###################################################################

###################################################################
# Konstrukcija posplošenega inverza preko Gramm-Schmidtove 
# ortogonalizacije
###################################################################

# norma vektorja
norma = function(x){
  sqrt(sum(x^2))
}

# GRAM-SCHMIDT -> Z = S * P
# S ima paroma pravokotne stolpce, nekateri so enaki 0
# P je zgornjetrikotna s pozitivnimi elementi na diagonali
GS = function(Z){
  n <- nrow(Z)
  d <- ncol(Z)
  S <- matrix(0, nrow = n, ncol = d)
  P <- matrix(0, nrow = d, ncol = d)
  # prvi stolpec S samo normiramo
  # drop = F uporabljamo zato, da se ohrani oblika (rezultat
  # bo še vedno matrika in ne vektor)
  S[, 1] = Z[, 1, drop = F] / norma(Z[, 1])
  P[1, 1] = norma(Z[, 1])
  
  for (i in 2:d){
    S[, i] = Z[, i, drop = F]
    # odštevamo pravokotne projekcije
    for (j in 1:(i-1)){
      koef = sum(Z[, i] * S[, j]) # skalarni produkt
      S[, i] = S[, i] - koef * S[, j]
      P[j, i] = koef
    }
    normaS = norma(S[,i])
    # če je pripadajoč stolpec enak 0, bo tudi norma enaka 0
    if (normaS < 1e-10){ 
      # na diagonalo P lahko postavimo katerokoli pozitivno število,
      # ker se množi z 0 (preko S)
      P[i, i] = 1
    }else{
      S[, i] = S[, i, drop = F] / normaS
      P[i, i] = normaS
    }
  }
  # vrnemo seznam z S in P
  list(S, P)
}


###################################################################
# Matrika Z, ki nima polnega ranga
# Datoteka primer1.csv na spletni učilnici
###################################################################

podatki = read.csv("primer1.csv")
Z = as.matrix(podatki[,3:6])
X = as.matrix(podatki[,2])

# Z = S * P
S = GS(Z)[[1]]
P = GS(Z)[[2]]

# Preverimo, če res dobimo Z, ko zmnožimo S in P
z = round(S %*% P, 10)
Z == z # vse OK :)

# še konstrukcija cenilke za beta
# lahko pa bi J dobili tudi preko norm stolpcev od S
J = t(S) %*% S
J = round(J, digits = 10)
# posplošeni inverz
inv_splosni = solve(P) %*% J %*% solve(t(P))

# cenilka
beta_hat = inv_splosni %*% t(Z) %*% X
beta_hat




