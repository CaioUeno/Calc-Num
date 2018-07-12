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

teste = function(A) {
    for (i in 1:nrow(A)) {
        if(A[i,i] == 0)
            return(0);
    }
    return(1);
}

Elim = function(A,B) {
  n = nrow(A);
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
    print(x);
    #print(y);
    e = max(dist)/max(abs(y));
    k = k + 1;
    #print(e);
  }
  cat("Erro da", k,"ªiteração:", e," \n", file = "SNL4.csv", append = TRUE);
  cat("\n", file = "SNL4.csv", append = TRUE);
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
  cat("Erro da", k,"ªiteração:", e," \n", file = "SNL4.csv", append = TRUE);
  cat("\n", file = "SNL4.csv", append = TRUE);
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
	if(inf == sup){ #https://gist.github.com/marcoscastro/60f8f82298212e267021
			cat(vetor,"\n", file="Permut_4.txt", sep=" ", append= T);
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
            return(i);
        }
    }
    return(0);
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
  W = c(0,0,0,0);
  cat("SNL\n", file = "SNL4.csv", append = FALSE); #limpa o arquivo
  while(E[1] > tol || E[2] > tol || E[3] > tol || E[4] > tol){
    c = c + 1;
    cat(c,"ª iteração:\n", file = "SNL4.csv", append = TRUE);
    A = Jacob(vetor_x); # aplica x na Jacobiana
    cat("Matriz Jacobiana aplicada em x:\n", file = "SNL4.csv", append = TRUE);
    write.table(A, file = "SNL4.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    B = matrix(c(-f1(vetor_x),-f2(vetor_x),-f3(vetor_x),-f4(vetor_x)),4,1);
    cat("Matriz B de - f's:\n", file = "SNL4.csv", append = TRUE);
    write.table(B, file = "SNL4.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    Flag = PermutaTCL(A,M_Permt); #Procura em todas as permutações possíveis
        if(Flag == 0){
            if(teste(A) == 1){
              vetor_a = Elim(A,B);
              cat("Resolvido por Eliminação.\n", file = "SNL4.csv", append = TRUE);
            } else {
              cat("Impossível resolver - (0 na diagonal principal!).\n", file = "SNL4.csv", append = TRUE);
              return(0);
            }
    }else{
        A = trocalinhas(A,M_Permt[,Flag]); #Permuta A, de forma a satisfazer o TCL
        B = trocalinhas(B,M_Permt[,Flag]); #Idem
        vetor_a = GaussJacobi(A,B,W,10^-6); # resolve A*vetor_a = B
        cat("Resolvido por Gauss-Jacobi.\n", file = "SNL4.csv", append = TRUE);
    }
    cat("Solução do Sistema Linear:\n", file = "SNL4.csv", append = TRUE);
    write.table(vetor_a, file = "SNL4.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    x_auxiliar = vetor_x; # armazena x(k)
    vetor_x = vetor_a + vetor_x; # x(k+1) = a + x(k)
    cat("Aproximação de x:\n", file = "SNL4.csv", append = TRUE);
    write.table(vetor_x, file = "SNL4.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
    E[1] = normaeuclidiana(vetor_x - x_auxiliar); #erro tipo 4
    E[2] = normaeuclidiana(F(vetor_x)); #erro tipo 5
    E[3] = max(abs(vetor_x - x_auxiliar)); #erro tipo 6
    E[4] = max(abs(F(vetor_x))); #erro tipo 7
    cat("Erros:\n", file = "SNL4.csv", append = TRUE);
    write.table(E, file = "SNL4.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
  }
  cat("Verificação da solução:\n", file = "SNL4.csv", append = TRUE);
  write.table(F(vetor_x), file = "SNL4.csv", sep = ", ", col.names = FALSE, row.names = FALSE, append = TRUE);
  cat("Tempo de execução:\n", file = "SNL4.csv", append = TRUE);
  }
normaeuclidiana = function(X) {
    s = 0;
    tam = length(X);
    for (i in 1:tam[1]) {
        s = s + X[i]^2;
    }
    return(sqrt(s));
}

Jacob = function(X) {
  MJ = cbind(dx1f(X),dx2f(X),dx3f(X),dx4f(X));
  print(MJ);
  return(MJ); # retorna a matriz Jacobiana do sistema
}
F = function(X) {
  Vec = c(0,0,0,0);
  Vec[1] = f1(X);
  Vec[2] = f2(X);
  Vec[3] = f3(X);
  Vec[4] = f4(X);
  return(Vec);
}

f2 = function(X){
    return( 16*X[1] + 2*X[2] - X[3]^2 + 4*X[4] - 16);
}
f1 = function(X){
    return( X[1]^2 + 19*X[2] - X[3] + cos(X[4]) - 42 );
}
f3 = function(X){
    return( 3*X[1] + 2*X[2] + 11*X[3]^4 - exp(X[4]) - 182 );
}
f4 = function(X){
    return( 5*X[1] + 3*X[2] - X[3] + 17*X[4] - 13 );
}
dx1f = function(X){
  der_f = matrix(c(0,0,0,0),4,1);
  der_f[2] = 16;
  der_f[1] = 2*X[1];
  der_f[3] = 3;
  der_f[4] = 5;
  return(der_f);
}

dx2f = function(X){
  der_f = matrix(c(0,0,0,0),4,1);
  der_f[2] = 2;
  der_f[1] = 19;
  der_f[3] = 2;
  der_f[4] = 3;
  return(der_f);
}

dx3f = function(X){
  der_f = matrix(c(0,0,0,0),4,1);
  der_f[2] = -2*X[3];
  der_f[1] = -1;
  der_f[3] = 11*4*X[3]^3;
  der_f[4] = -1;
  return(der_f);
}

dx4f = function(X){
  der_f = matrix(c(0,0,0,0),4,1);
  der_f[2] = 4;
  der_f[1] = -sin(X[4]);
  der_f[3] = -exp(X[4]);
  der_f[4] = 17;
  return(der_f);
}
Q = c(1,2,-2,0);#solução
X = c(5,8,-4,2);
X = c(10,10,30,10);
X = c(50,80,-100,-45);
vetu = c(1,2,3,4);
#permuta(A,vetu, 1,4);
M_Permt = matrix(scan("Permut_4.txt"),4,factorial(4));
SistemaNL(X,10^-12);

print(proc.time());
