TCL = function(A){
  n = nrow(A);
  for (i in 1:n) {
    s = 0;
    for (j in 1:n) {
      if(i!=j){
        s = s + abs(A[i,j]);
      }
    }
    s = s/abs(A[i,i]);
    if(s >= 1){
      return (0);
    }
  }
  return (1);
}

Sassenfeld = function(A){
  n = nrow(A);
  Beta = c(0);
  for(i in 1:n){
    s = 0;
    for (j in 1:n) {
      if(i!=j){
        if(j<i){
          s = s + abs(A[i,j])*Beta[j];
        } else {
          s = s + abs(A[i,j]);
        }
      }
    }
    Beta[i] = s/abs(A[i,i]);
    if(is.infinite(Beta[i]) || is.nan(Beta[i])){
        return(0);
    }
  }
  if(max(Beta) < 1){
    print(Beta);
    return(1);
  } else {
    return(0);
  }
}

teste = function(A) { #checa se não há 0's na diagonal principal
    for(i in nrow(A)){
        if(A[i,i]==0){
            return(0);
        }
    }
    return(1);
}

Elim = function(A,B) {
  n = nrow(A);
  X = c(0);
  for(k in 1:(n-1)) {
    for(i in (k+1):n) {
      m = A[i,k]/A[k,k];
      A[i,k] = 0;
      for(j in (k+1):n) {
        A[i,j] = A[i,j] - m*A[k,j];
      }
      B[i] = B[i] - m*B[k];
    }
  }

  X[n] = B[n] / A[n, n]

  for (i in (n - 1):1) {
    soma = 0

    for (j in (i + 1):n) {
      soma = soma + A[i, j] * X[j]

    }
    X[i] = (B[i] - soma) / A[i, i]
  }

  return(X);
}

GaussJacobi = function(A, B, C, tol){
  n = nrow(A);
  x = C;
  dist = c(0);
  e = 1;
  k = 1;
  #print(A);
  #print(B);
  while(e > tol){
    y = c(0);
    for(i in 1:n){
      s = 0;
      for(j in 1:n){
        if(i!=j){
          s = s + A[i,j]*x[j];
        }
      }
      y[i] = (B[i] - s)/A[i,i];
      dist[i] = abs(y[i] - x[i]);
    }
    x = y;
    #print(x);
    #print(y);
    e = max(dist)/max(abs(y));
    k = k + 1;
    #print(e);
  }
  cat("Erro da", k," iteracao:", e," \n", file = "SNL6.csv", append = TRUE);
  cat("\n", file = "SNL6.csv", append = TRUE);
  #print(x);
  #print(k);
  return(x);
}

GaussSeidel = function(A, B, C, tol){
  n = nrow(A);
  x = C;
  dist = c(0);
  e = 1;
  k = 1;
  while(e > tol){
    y = c(0);
    for(i in 1:n){
      s = 0;
      for(j in 1:n){
        if(i!=j){
          s = s + A[i,j]*x[j];
        }
      }
      y[i] = (B[i] - s)/A[i,i];
      dist[i] = abs(y[i] - x[i]);
      x[i] = y[i];
    }
    #print(x);
    #print(y);
    e = max(dist)/max(abs(y));
    print(e);
    k = k + 1;
    #print(e);
  }
  cat("Erro da", k,"  iteracao:", e," \n", file = "SNL6.csv", append = TRUE);
  cat("\n", file = "SNL6.csv", append = TRUE);
  #print(x);
  #print(k);
  return(x);
}

troca = function(vetor,i,j){
  aux = vetor[i];
  vetor[i] = vetor[j];
  vetor[j] = aux;
  return(vetor);
}

permuta = function(A,vetor,inf,sup){ #Cria arquivo com todas as permutacoes de n elementos
  if(inf == sup){
    cat(vetor,"\n", file="Permut_6.txt", sep=" ", append= T);
  }
  else{
    for(i in inf:sup){
      vetor = troca(vetor, inf, i);
      permuta(A,vetor, inf + 1, sup);
      vetor = troca(vetor, inf, i); # backtracking
    }
  }
}

trocalinhas = function(A,ordem){ #Troca as linhas de uma matriz numa ordem escolhida
  B = A;
  for (i in 1:nrow(A)) {
    A[i,] = B[ordem[i],];
  }
  return(A);
}

PermutaTCL = function(A,MP){ # testa o TCL em todas as permutacos possiveis de linhas
  for (i in 1:ncol(MP)) {
    H = trocalinhas(A,MP[,i]);
    if(TCL(H)==1){
      return(i); #retorna a coluna em que tem a permutacao que satisfaz o TCL
    }
  }
  return(0); # caso contrario retorna 0
}

PermutaSassenfeld = function(A,MP){ # testa o Sassenfeld em todas as permutacoes possiveis de linhas
  for (i in 1:ncol(MP)) {
    H = trocalinhas(A,MP[,i]);
    if(Sassenfeld(H)==1){
      return(i); #retorna a coluna em que tem a permutacao que satisfaz o Sassenfeld
    }
  }
  return(0); # caso contrario retorna 0
}

SistemaNL = function(vetor_x, tol) {
  vetor_a = c(0);
  E = c(1,1,1,1); # vetor com os 4 tipos de erro de SNL
  c = 0;
  W = c(0,0,0,0,0,0); # chute inicial para Jacobi/Seidel - esperado para se chegar na solucao
  cat("SNL\n", file = "SNL6.csv", append = FALSE);
  while(E[1] > tol || E[2] > tol || E[3] > tol || E[4] > tol){
    c = c + 1;
    cat(c,"ª iteracao:\n", file = "SNL6.csv", append = TRUE);
    A = Jacob(vetor_x); # aplica x na Jacobiana
    cat("Matriz Jacobiana aplicada em x:\n", file = "SNL6.csv", append = TRUE);
    write.table(A, file = "SNL6.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    B = matrix(c(-f1(vetor_x),-f2(vetor_x),-f3(vetor_x),-f4(vetor_x),-f5(vetor_x),-f6(vetor_x)),6,1);
    print(A);
    cat("Matriz B de - f's:\n", file = "SNL6.csv", append = TRUE);
    write.table(B, file = "SNL6.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    #Flag = 0;
    Flag = PermutaSassenfeld(A,M_Permt); #Procura em todas as permutacoes possi?veis
    if(Flag == 0){
      if(teste(A)==0){
          cat("Não é possível resolver por Eliminação.\n", file = "SNL6.csv", append = TRUE);
          return(3);
      } else{
          vetor_a = Elim(A,B);
          cat("Resolvido por Eliminação.\n", file = "SNL6.csv", append = TRUE);
      }
    }else{
      A = trocalinhas(A,M_Permt[,Flag]); #Permuta A, de forma a satisfazer o TCL
      B = trocalinhas(B,M_Permt[,Flag]); #Idem
      vetor_a = GaussSeidel(A,B,W,10^-6); # resolve A*vetor_a = B
      cat("Resolvido por Gauss-Seidel.\n", file = "SNL6.csv", append = TRUE);
    }
    cat("Solucao do Sistema Linear:\n", file = "SNL6.csv", append = TRUE);
    write.table(vetor_a, file = "SNL6.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    x_auxiliar = vetor_x; # armazena x(k)
    vetor_x = vetor_a + vetor_x; # x(k+1) = a + x(k)
    cat("Aproximacao de x:\n", file = "SNL6.csv", append = TRUE);
    write.table(vetor_x, file = "SNL6.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    E[1] = normaeuclidiana(vetor_x - x_auxiliar); #erro tipo 4
    E[2] = normaeuclidiana(F(vetor_x)); #erro tipo 5
    E[3] = max(abs(vetor_x - x_auxiliar)); #erro tipo 6
    E[4] = max(abs(F(vetor_x))); #erro tipo 7
    cat("Erros:\n", file = "SNL6.csv", append = TRUE);
    write.table(E, file = "SNL6.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
  }
  cat("Verificacao da solucao:\n", file = "SNL6.csv", append = TRUE);
  write.table(F(vetor_x), file = "SNL6.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);

}
normaeuclidiana = function(X) {
  s = 0;
  tam = length(X);
  #print(tam[1]);
  for (i in 1:tam[1]) {
    s = s + X[i]^2;
  }
  return(sqrt(s));
}

Jacob = function(X) { #concatena os vetores - em coluna - para construir a Jacobiana
  MJ = cbind(dx1f(X),dx2f(X),dx3f(X),dx4f(X),dx5f(X),dx6f(X));
  print(MJ);
  return(MJ); # retorna a matriz Jacobiana do sistema
}
#Aplica x em todas as funÃ§Ãµes, para fins de teste
F = function(X) {
  Vec = c(0,0,0,0,0);
  Vec[1] = f1(X);
  Vec[2] = f2(X);
  Vec[3] = f3(X);
  Vec[4] = f4(X);
  Vec[5] = f5(X);
  Vec[6] = f6(X);
  return(Vec);
}

f1 = function(X){
  return(159*X[1] + 2*X[2]^2 + 5*X[3] - 2*X[4]/5 - 3*X[5]^3 - 3*X[6] -183);
}
f2 = function(X){
  return(X[1]^3 - 4*X[2]^4 - cos(X[3]) + 2*X[4]^2 + X[5]^3  - X[6] + 24);
}
f3 = function(X){
  return(2*X[1] + 3*X[2] - 27*X[3] - 1/4*X[4]^2 - X[5]^2 + 0.01*X[6]^3 + 2.17);
}
f4 = function(X){
  return(-1*X[1]^2 - X[2] + 1.4*sin(X[3]) - 5*(X[4]^3) + (X[5])/4 + 1.1*X[6]^2 + 624.10);
}
f5 = function(X) {
  return(2*X[1] - 2*X[2] + 0*X[3] - 5*X[4] + 21*X[5] - 2*X[6] + 73);
}
f6 = function(X) {
  return(-8*X[1] - X[2]/2 - 0.01*exp(X[3]) + 2*X[4]/5 + X[5] + 9*X[6]^5 -278.99);
}
#As funcoes seguintes retornam colunas da matriz Jacobiana, para diminuir a quantidade de parametros
dx1f = function(X){ #retorna um vetor com as derivadas em relacao a x1 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0),6,1);
  der_f[1] = 159;
  der_f[2] = 3*X[1]^2;
  der_f[3] = 2;
  der_f[4] = -2*X[1];
  der_f[5] = 2;
  der_f[6] = -8;
  return(der_f);
}
dx2f = function(X){ #retorna um vetor com as derivadas em relacao a x2 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0),6,1);
  der_f[1] = 4*X[2];
  der_f[2] = -16*X[2]^3;
  der_f[3] = 3;
  der_f[4] = -1;
  der_f[5] = -2;
  der_f[6] = -1/2;
  return(der_f);
}
dx3f = function(X){ #retorna um vetor com as derivadas em relacao a x3 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0),6,1);
  der_f[1] = 5;
  der_f[2] = sin(X[3]);
  der_f[3] = -27;
  der_f[4] = 1.4*cos(X[3]);
  der_f[5] = 0;
  der_f[6] = -0.01*exp(X[3]);
  return(der_f);
}
dx4f = function(X){ #retorna um vetor com as derivadas em relacao a x4 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0),6,1);
  der_f[1] = -2/5;
  der_f[2] = 4*X[4];
  der_f[3] = -0.5*X[4];
  der_f[4] = -15*X[4]^2;
  der_f[5] = 5;
  der_f[6] = 2/5;
  return(der_f);
}
dx5f = function(X){ #retorna um vetor com as derivadas em relacao a x5 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0),6,1);
  der_f[1] = -9*X[5]^2;
  der_f[2] = 3*X[5]^2;
  der_f[3] = -2*X[5];
  der_f[4] = 1/4;
  der_f[5] = 21;
  der_f[6] = 1;
  return(der_f);
}
dx6f = function(X){ #retorna um vetor com as derivadas em relacao a x6 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0),6,1);
  der_f[1] = -3;
  der_f[2] = -1;
  der_f[3] = 0.03*X[6]^2;
  der_f[4] = 2.2*X[6];
  der_f[5] = -2;
  der_f[6] = 45*X[6]^4;
  return(der_f);
}
C = c(1,2,0,5,-2,2);#solução
C = c(8,5,3,11,-8,5);#ok
C = c(20,21,20,20,22,22);#ok
C = c(1003,2,0,5,-2,2);#nao converge pra solução exata
#Jacob(S);
F(C);
#vetu = c(1,2,3,4,5,6);
#permuta(A,vetu, 1,6);
M_Permt = matrix(scan("Permut_6.txt"),6,factorial(6));
# matriz em que cada coluna possui uma permutacao de 1 a n, no total de n! colunas
SistemaNL(C,10^-12);
print(proc.time());
