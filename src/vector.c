#include "vector.h"
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <errno.h>
#include <stddef.h>
#include <math.h>

Vector* vector_create_empty(unsigned int dim) {
  if (dim == 0) {
    errno = EINVALIDDIM;
    return NULL;
  }
  
  Vector* v = (Vector*)malloc(sizeof(Vector));
  if (v == NULL) {
    errno = EVECTORALLOCATION;
    return NULL;
  }

  v->dim = dim;
  v->x = (double*)calloc(dim, sizeof(double));
  if (v->x == NULL) {
    errno = EVECTORALLOCATION;
    free(v);
    return NULL;
  }

  return v;
}

Vector* vector_create(unsigned int dim, double coords[]) {
  if (dim == 0) {
    errno = EINVALIDDIM;
    return NULL;
  }

  Vector* v = (Vector*)malloc(sizeof(Vector));
  if (v == NULL) {
    errno = EVECTORALLOCATION;
    return NULL;
  }

  v->dim = dim;
  v->x = (double*)calloc(dim, sizeof(double));
  if (v->x == NULL) {
    errno = EVECTORALLOCATION;
    free(v);
    return NULL;
  }

  memcpy(v->x, coords, dim * sizeof(double));

  return v;
}

void vector_destroy(Vector* v) {
  if (v != NULL) {
    if (v->x != NULL) {
      free(v->x);
    }
    free(v);
  }
}

Vector* vector_add(const Vector* a, const Vector* b) {
  if (a == NULL || b == NULL) {
    errno = EVECTORNULL;
    return NULL;
  }
  
  if (a->dim != b->dim) {
    errno = EINVALIDDIM;
    return NULL;
  }

  Vector* v = (Vector*)malloc(sizeof(Vector));
  if (v == NULL) {
    errno = EVECTORALLOCATION;
    return NULL;
  }

  v->dim = a->dim;
  v->x = (double*)calloc(a->dim, sizeof(double));
  if (v->x == NULL) {
    errno = EVECTORALLOCATION;
    free(v);
    return NULL;
  }

  for (unsigned int i = 0; i < a->dim; i++) {
    v->x[i] = a->x[i] + b->x[i];
  }
  
  return v;
}

Vector* vector_sub(const Vector* a, const Vector* b) {
  if (a == NULL || b == NULL) {
    errno = EVECTORNULL;
    return NULL;
  }
 
  if (a->dim != b->dim) {
    errno = EINVALIDDIM;
    return NULL;
  }
  
  Vector* v = (Vector*)malloc(sizeof(Vector));
  if (v == NULL) {
    errno = EVECTORALLOCATION;
    return NULL;
  }
  
  v->dim = a->dim;
  v->x = (double*)calloc(a->dim, sizeof(double));
  if (v->x == NULL) {
    errno = EVECTORALLOCATION;
    free(v);
    return NULL;
  }

  for (unsigned int i = 0; i < a->dim; i++) {
    v->x[i] = a->x[i] - b->x[i];
  }
  
  return v;
}

Vector* vector_mul(double scalar, const Vector* a) {
  if (a == NULL) {
    errno = EVECTORNULL;
    return NULL;
  }
  
  Vector* v = (Vector*)malloc(sizeof(Vector));
  if (v == NULL) {
    errno = EVECTORALLOCATION;
    return NULL;
  }

  v->dim = a->dim;
  v->x = (double*)calloc(a->dim, sizeof(double));
  if (v->x == NULL) {
    errno = EVECTORALLOCATION;
    free(v);
    return NULL;
  }

  for (unsigned int i = 0; i < a->dim; i++) {
    v->x[i] = scalar * a->x[i];
  }
  
  return v;
}

double vector_inner_product(const Vector* a, const Vector* b) {
  if (a == NULL || b == NULL) {
    errno = EVECTORNULL;
    return NAN;
  }

  if (a->dim != b->dim) {
    errno = EINVALIDDIM;
    return NAN;
  }

  double p = 0.00;
  for (unsigned int i = 0; i< a->dim; i++) {
    p += a->x[i]* b->x[i];
  }

  return p;
}

double vector_norm_p(const Vector* a, unsigned int p) {
  if (a == NULL) {
    errno = EVECTORNULL;
    return NAN;
  }

  if (p == 0) {
    errno = EINVALIDP;
    return NAN;
  }

  double sum = 0.00;
  for (unsigned int i = 0; i < a->dim; i++) {
    sum += pow(fabs(a->x[i]), (double)p);
  }
  double rad = 1.0 / (double)p;

  return pow(sum, rad);
}

double vector_norm_inf(const Vector* a) {
  if (a == NULL) {
    errno = EVECTORNULL;
    return NAN;
  }

  double max = fabs(a->x[0]);
  for (unsigned int i = 1; i < a->dim; i++) {
    double modulo = fabs(a->x[i]);
    if (modulo > max) {
      max = modulo;
    }
  }

  return max;
}

double vector_norm(const Vector* a) {
  return vector_norm_p(a, 2);
}

Vector* vector_product(const Vector* a, const Vector* b) {
//Algoritmo "Produto_Vetorial_Dimensional_N"
//// Descrição   : Calcula o vetor resultante (ortogonal) a partir de (n-1) vetores de dimensão n
//// Vetores de entrada: M[1..n-1, 1..n]
//// Vetor de saída: VetorResultante[1..n]
//
//Var
//   n, i, j: inteiro
//   M: vetor [1..20, 1..20] de real // Matriz para armazenar os n-1 vetores
//   VetorResultante: vetor [1..20] de real
//   SubMatriz: vetor [1..20, 1..20] de real
//   detSub: real
//
//// Função auxiliar fictícia para calcular o determinante de uma matriz quadrada de ordem (n-1)
//// Em um sistema real, esta função implementaria o método de Gauss ou Teorema de Laplace recursivo.
//Funcao CalcularDeterminante(Mat: vetor [1..20, 1..20] de real; ordem: inteiro): real
//Var
//   // Variáveis de controle do determinante
//Inicio
//   // O cálculo real do determinante de ordem (n-1) entraria aqui.
//   // Para n=3 (submatriz 2x2): retorna (Mat[1,1]*Mat[2,2] - Mat[1,2]*Mat[2,1])
// FimFuncao
//
//Inicio
//   Escreva("Digite a dimensão N do espaço (ex: 3, 4, 5): ")
//   Leia(n)
//   
//   Se n < 2 Entao
//      Escreval("A dimensão deve ser de pelo menos 2.")
//      FimAlgoritmo
//   FimSe
//
//   Escreval("Para o produto vetorial em ", n, " dimensões, você precisa fornecer ", n-1, " vetores.")
//   
//   // Leitura dos n-1 vetores
//   Para i de 1 ate n-1 faca
//      Escreval("--- Vetor ", i, " ---")
//      Para j de 1 ate n faca
//         Escreva("Componente ", j, ": ")
//         Leia(M[i, j])
//      FimPara
//   FimPara
//
//   // Cálculo do "produto vetorial" via Cofatores da primeira linha imaginária
//   Para j de 1 ate n faca
//      // Monta a submatriz eliminando a coluna j
//      ConstruirSubMatriz(M, SubMatriz, n, j)
//      
//      // Calcula o determinante da submatriz de ordem n-1
//      detSub <- CalcularDeterminante(SubMatriz, n-1)
//      
//      // Aplica o sinal do cofator: (-1)^(1 + j) -> se j for par, inverte o sinal
//      Se (1 + j) % 2 = 0 Entao
//         VetorResultante[j] <- detSub
//      Senao
//         VetorResultante[j] <- -detSub
//      FimSe
//   FimPara
//
//   // Exibe o vetor ortogonal resultante
//   Escreval("---------------------------------------")
//   Escreva("Vetor Resultante (Produto Vetorial): [")
//   Para j de 1 ate n faca
//      Escreva(VetorResultante[j])
//      Se j < n Entao
//         Escreva(", ")
//      FimSe
//   FimPara
//   Escreval("]")
//
//FimAlgoritmo
//
//// Procedimento para isolar a submatriz eliminando a coluna 'colunaExcluida'
//Procedimento ConstruirSubMatriz(MatOrigem: vetor [1..20, 1..20] de real; Var MatDestino: vetor [1..20, 1..20] de real; dimn, colunaExcluida: inteiro)
//Var
//   i, j, colAtual: inteiro
//Inicio
//   Para i de 1 ate dimn-1 faca
//      colAtual <- 1
//      Para j de 1 ate dimn faca
//         Se j <> colunaExcluida Entao
//            MatDestino[i, colAtual] <- MatOrigem[i, j]
//            colAtual <- colAtual + 1
//         FimSe
//      FimPara
//   FimPara
//FimProcedimento

}

void vector_print(const Vector* v, const char* vector_label) {
  if (v == NULL) {
    errno = EVECTORNULL;
    fprintf(stderr, "==>ERRO: EVECTORNULL");
    return;
  }

  if (vector_label != NULL) {
    printf("%s = ", vector_label);
  }

  printf("(");
  for (unsigned int i = 0; i < v->dim; i++) {
    if (i + 1 < v->dim) {
      printf("%.2f, ", v->x[i]);
    } else {
      printf("%.2f", v->x[i]);
    }
  }
  printf(")");
}

void vector_println(const Vector* v, const char* vector_label) {
  vector_print(v, vector_label);
  printf("\n");
}

