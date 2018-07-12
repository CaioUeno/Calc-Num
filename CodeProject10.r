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
  cat("Erro da", k,"ªiteração:", e," \n", file = "SNL10.csv", append = TRUE);
  cat("\n", file = "SNL10.csv", append = TRUE);
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
    k = k + 1;
    #print(e);
  }
  cat("Erro da", k,"ªiteração:", e," \n", file = "SNL10.csv", append = TRUE);
  cat("\n", file = "SNL10.csv", append = TRUE);
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

permuta = function(A,vetor,inf,sup){ #Cria arquivo com todas as permutações de n elementos
  if(inf == sup){
    cat(vetor,"\n", file="Permut_10.txt", sep=" ", append= T);
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

PermutaTCL = function(A,MP){ # testa o TCL em todas as permutações possíveis de linhas
  for (i in 1:ncol(MP)) {
    H = trocalinhas(A,MP[,i]);
    if(TCL(H)==1){
      return(i); #retorna a coluna em que tem a permutação que satisfaz o TCL
    }
  }
  return(0); # caso contrário retorna 0
}

PermutaSassenfeld = function(A,MP){ # testa o Sassenfeld em todas as permutações possíveis de linhas
  for (i in 1:ncol(MP)) {
    H = trocalinhas(A,MP[,i]);
    if(Sassenfeld(H)==1){
      return(i); #retorna a coluna em que tem a permutação que satisfaz o Sassenfeld
    }
  }
  return(0); # caso contrário retorna 0
}

SistemaNL = function(vetor_x, tol) {
  vetor_a = c(0);
  E = c(1,1,1,1); # vetor com os 4 tipos de erro de SNL
  c = 0;
  W = c(0,0,0,0,0,0,0,0,0,0); # chute inicial para Jacobi/Seidel - esperado para se chegar na solução
  cat("SNL\n", file = "SNL10.csv", append = FALSE);
  while(E[1] > tol || E[2] > tol || E[3] > tol || E[4] > tol){
    c = c + 1;
    cat(c,"ª iteração:\n", file = "SNL10.csv", append = TRUE);
    A = Jacob(vetor_x); # aplica x na Jacobiana
    cat("Matriz Jacobiana aplicada em x:\n", file = "SNL10.csv", append = TRUE);
    write.table(A, file = "SNL10.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    B = matrix(c(-f1(vetor_x),-f2(vetor_x),-f3(vetor_x),-f4(vetor_x),-f5(vetor_x),
                 -f6(vetor_x),-f7(vetor_x),-f8(vetor_x),-f9(vetor_x),-f10(vetor_x)),10,1);
    cat("Matriz B de - f's:\n", file = "SNL10.csv", append = TRUE);
    write.table(B, file = "SNL10.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    #Flag = PermutaTCL(A,M_Permt); #Procura em todas as permutações possíveis
    #if(Flag == 0){
    if(teste(A)==0){
        cat("Não é possível resolver por Eliminação.\n", file = "SNL10.csv", append = TRUE);
        return(3);
    } else{
        vetor_a = Elim(A,B);
        cat("Resolvido por Eliminação.\n", file = "SNL10.csv", append = TRUE);
    }
    #}else{
    #  A = trocalinhas(A,M_Permt[,Flag]); #Permuta A, de forma a satisfazer o TCL
    #  B = trocalinhas(B,M_Permt[,Flag]); #Idem
    #  vetor_a = GaussJacobi(A,B,W,10^-12); # resolve A*vetor_a = B
    #  cat("Resolvido por Gauss-Jacobi.\n", file = "SNL10.csv", append = TRUE);
    #}
    cat("Solução do Sistema Linear:\n", file = "SNL10.csv", append = TRUE);
    write.table(vetor_a, file = "SNL10.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    x_auxiliar = vetor_x; # armazena x(k)
    vetor_x = vetor_a + vetor_x; # x(k+1) = a + x(k)
    cat("Aproximação de x:\n", file = "SNL10.csv", append = TRUE);
    write.table(vetor_x, file = "SNL10.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    E[1] = normaeuclidiana(vetor_x - x_auxiliar); #erro tipo 4
    E[2] = normaeuclidiana(F(vetor_x)); #erro tipo 5
    E[3] = max(abs(vetor_x - x_auxiliar)); #erro tipo 6
    E[4] = max(abs(F(vetor_x))); #erro tipo 7
    cat("Erros:\n", file = "SNL10.csv", append = TRUE);
    write.table(E, file = "SNL10.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
  }
  cat("Verificação da solução:\n", file = "SNL10.csv", append = TRUE);
  write.table(F(vetor_x), file = "SNL10.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);

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
  MJ = cbind(dx1f(X),dx2f(X),dx3f(X),dx4f(X),dx5f(X),
             dx6f(X),dx7f(X),dx8f(X),dx9f(X),dx10f(X));
  print(MJ);
  return(MJ); # retorna a matriz Jacobiana do sistema
}
#Aplica x em todas as funções, para fins de teste
F = function(X) {
  Vec = c(0,0,0,0,0,0,0,0,0,0);
  Vec[1] = f1(X);
  Vec[2] = f2(X);
  Vec[3] = f3(X);
  Vec[4] = f4(X);
  Vec[5] = f5(X);
  Vec[6] = f6(X);
  Vec[7] = f7(X);
  Vec[8] = f8(X);
  Vec[9] = f9(X);
  Vec[10] = f10(X);
  return(Vec);
}

f1 = function(X){
  return(X[1] + X[2] + X[3] + X[4] + X[5] + X[6] + X[7] + X[8] + X[9] + X[10] + 2);
}
f2 = function(X){
  return(X[1] + X[2]^2 + X[3]^3 + X[4]^4 + X[5]^5 + X[6]^4 + X[7]^3 + X[8]^2 + X[9] + 5*X[10] - 56);
}
f3 = function(X){
  return(10*X[1] - X[2] + exp(X[3]) + 3*X[4]^2 - X[5] + 2*X[6]^3 - X[7]^2 + 2*X[10] + 48);
}
f4 = function(X){
  return(2*X[2] + cos(X[3]) + sin(X[4]) - 3*X[5] + (2*X[6])^2 - 12*X[7] + exp(X[8]) + 3*X[9] + 45);
}
f5 = function(X) {
  return(-2*X[1] + 3*(X[2]^3) - exp(X[4]) + X[5]^2 + 3*X[6] - X[7]^2 - X[9]^3 + X[10] + 70);
}
f6 = function(X) {
  return(-X[1]^2 + X[2]^4 + sin(X[3]) - X[4] + 10*X[5] - 2*X[6] - X[7]^2 - 2^(X[8]) - X[9]^3 + 2*X[10] - 76);
}
f7 = function(X) {
  return(-3*X[1] + cos(X[2] + 3) - X[4] + X[5]^2 - X[6] - X[7]^3 - X[9] - X[10] + 12);
}
f8 = function(X) {
  return(-X[1] - X[2]^4 + 3*cos(X[3]) - X[4]^3 + X[5]^5 - X[6] + 18*X[7] - 0.7*exp(X[8]) + X[9] + 2*X[10] + 9.7);
}
f9 = function(X) {
  return(2*cos(X[3]) - 0.1*(X[4]) + 2*X[5] - X[6]^3 - 5*exp(X[8]) );
}
f10 = function(X) {
  return(X[5] - cos(X[8]) + 40*X[10] - 200);
}

#As funções seguintes retornam colunas da matriz Jacobiana, para diminuir a quantidade de parâmetros
dx1f = function(X){ #retorna um vetor com as derivadas em relação a x1 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0,0,0,0,0),10,1);
  der_f[1] = 1;
  der_f[2] = 1;
  der_f[3] = 10;
  der_f[4] = 0
  der_f[5] = -2;
  der_f[6] = -2*X[1];
  der_f[7] = -3;
  der_f[8] = -1;
  der_f[9] = 0;
  der_f[10] = 0;
  return(der_f);
}
dx2f = function(X){ #retorna um vetor com as derivadas em relação a x2 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0,0,0,0,0),10,1);
  der_f[1] = 1;
  der_f[2] = 2*X[2];
  der_f[3] = -1;
  der_f[4] = 2;
  der_f[5] = 9*X[2]^2;
  der_f[6] = 4*X[2]^3;
  der_f[7] = -sin(X[2] + 3);
  der_f[8] = -4*X[2]^3;
  der_f[9] = 0;
  der_f[10] = 0;
  return(der_f);
}
dx3f = function(X){ #retorna um vetor com as derivadas em relação a x3 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0,0,0,0,0),10,1);
  der_f[1] = 1;
  der_f[2] = 3*X[2]^2;
  der_f[3] = exp(X[3]);
  der_f[4] = -sin(X[3]);
  der_f[5] = 0;
  der_f[6] = cos(X[3]);
  der_f[7] = 0;
  der_f[8] = -3*sin(X[3]);
  der_f[9] = -2*sin(X[3]);
  der_f[10] = 0;
  return(der_f);
}
dx4f = function(X){ #retorna um vetor com as derivadas em relação a x4 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0,0,0,0,0),10,1);
  der_f[1] = 1;
  der_f[2] = 4*X[4]^3;
  der_f[3] = 6*X[4];
  der_f[4] = cos(X[4]);
  der_f[5] = -exp(X[4]);
  der_f[6] = -1;
  der_f[7] = -1;
  der_f[8] = -3*X[4]^2;
  der_f[9] = -0.1;
  der_f[10] = 0;
  return(der_f);
}
dx5f = function(X){ #retorna um vetor com as derivadas em relação a x5 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0,0,0,0,0),10,1);
  der_f[1] = 1;
  der_f[2] = 5*X[5]^4;
  der_f[3] = -1;
  der_f[4] = -3;
  der_f[5] = 2*X[5];
  der_f[6] = 10;
  der_f[7] = 2*X[5];
  der_f[8] = 5*X[5]^4;
  der_f[9] = 2;
  der_f[10] = 1;
  return(der_f);
}
dx6f = function(X){ #retorna um vetor com as derivadas em relação a x6 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0,0,0,0,0),10,1);
  der_f[1] = 1;
  der_f[2] = 4*X[6]^3;
  der_f[3] = 6*X[6]^2;
  der_f[4] = 4*X[6];
  der_f[5] = 3;
  der_f[6] = -2;
  der_f[7] = -1;
  der_f[8] = -1;
  der_f[9] = -3*X[6]^2;
  der_f[10] = 0;
  return(der_f);
}
dx7f = function(X){ #retorna um vetor com as derivadas em relação a x7 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0,0,0,0,0),10,1);
  der_f[1] = 1;
  der_f[2] = 3*X[7]^2;
  der_f[3] = -2*X[7];
  der_f[4] = -12;
  der_f[5] = -2*X[7];
  der_f[6] = -2*X[7];
  der_f[7] = -3*X[7]^2;
  der_f[8] = 18;
  der_f[9] = 0;
  der_f[10] = 0;
  return(der_f);
}
dx8f = function(X){ #retorna um vetor com as derivadas em relação a x8 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0,0,0,0,0),10,1);
  der_f[1] = 1;
  der_f[2] = 2*X[8];
  der_f[3] = 0;
  der_f[4] = exp(X[8]);
  der_f[5] = 0;
  der_f[6] = (-1)*log(2)*(2^X[8]);
  der_f[7] = 0;
  der_f[8] = -0.7*exp(X[8]);
  der_f[9] = -5*exp(X[8]);
  der_f[10] = sin(X[8]);
  return(der_f);
}
dx9f = function(X){ #retorna um vetor com as derivadas em relação a x9 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0,0,0,0,0),10,1);
  der_f[1] = 1;
  der_f[2] = 1;
  der_f[3] = 0;
  der_f[4] = 3;
  der_f[5] = -3*X[9]^2;
  der_f[6] = -3*X[9]^2;
  der_f[7] = -1;
  der_f[8] = 1;
  der_f[9] = 0;
  der_f[10] = 0;
  return(der_f);
}
dx10f = function(X){ #retorna um vetor com as derivadas em relação a x10 de todas as f's
  der_f = matrix(c(0,0,0,0,0,0,0,0,0,0),10,1);
  der_f[1] = 1;
  der_f[2] = 5;
  der_f[3] = 2;
  der_f[4] = 0;
  der_f[5] = 1;
  der_f[6] = 2;
  der_f[7] = -1;
  der_f[8] = 2;
  der_f[9] = 0;
  der_f[10] = 40;
  return(der_f);
}
C = c(-5,-3,0,0,1,-1,3,0,-2,5);
C = c(-5.00001,-3.00002,0.00003,0.001,1.0003,-1.00000002,3.00007,0.0000002,-2,5);
C =c(-5.9,-3.8,0.9,0,1,-1,3,0.9,-2.7,5.1);
#C =c(-7,-4,4,4,3,-3.5,3.7,1.3,-5,9); # erro
#Jacob(S);
F(C);
#vetu = c(1,2,3,4,5,6,7,8,9,10);
#permuta(A,vetu, 1,10);
#M_Permt = matrix(scan("Permut_10.txt"),10,factorial(10));
# matriz em que cada coluna possui uma permutação de 1 a n, no total de n! colunas
SistemaNL(C,10^-12);
