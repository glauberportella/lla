#include "lu.h"
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <math.h>

/**
 * Cria a matriz de permutacao P a partir do vetor de pivots.
 *
 * O vetor p indica que a linha i de P*A corresponde a linha p[i] de A,
 * ou seja, P[i][p[i]] = 1 e as demais entradas da linha i sao zero.
 */
static Matrix* lu_create_permutation(const unsigned int n, const int p[]) {
  Matrix* P = matrix_create_empty(n, n);
  if (P == NULL) {
    return NULL;
  }

  for (unsigned int i = 0; i < n; i++) {
    matrix_set(P, i, (unsigned int)p[i], 1.0);
  }

  return P;
}

LU* lu_create_decomposition(const Matrix* a) {
  if (a == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  if (a->nlins != a->ncols) {
    errno = EMATRIXINVALIDDIM;
    return NULL;
  }

  if (a->nlins == 0) {
    errno = EMATRIXINVALIDDIM;
    return NULL;
  }

  Matrix* L = matrix_create_empty(a->nlins, a->ncols);
  if (L == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  Matrix* U = matrix_create_empty(a->nlins, a->ncols);
  if (U == NULL) {
    matrix_destroy(L);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  // Vetor de permutação que registra as trocas de linha
  int* P = (int*)calloc(a->nlins, sizeof(int));
  if (P == NULL) {
    matrix_destroy(U);
    matrix_destroy(L);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  unsigned int i, j, k, linha_max;
  double max_val, aux, fator;

  // 1. Inicializar matrizes e o vetor de permutação P
  for (i = 0; i < a->nlins; i++) {
    P[i] = (int)i;
    // Diagonal de L sempre com 1.0
    matrix_set(L, i, i, 1.0);
  }

  // Copiar a matriz A para U, pois faremos o escalonamento em U
  for (i = 0; i < a->nlins; i++) {
    for (j = 0; j < a->ncols; j++) {
      matrix_set(U, i, j, matrix_get(a, i, j));
    }
  }

  // 2. Processo de Decomposição LU com Pivoteamento Parcial
  for (i = 0; i < a->nlins; i++) {
    // --- PASSO DE PIVOTEAMENTO PARCIAL ---
    linha_max = i;
    max_val = fabs(matrix_get(U, i, i));

    // Procura pelo maior elemento em módulo na coluna i, da linha i para baixo
    for (k = i + 1; k < a->nlins; k++) {
      if (fabs(matrix_get(U, k, i)) > max_val) {
        max_val = fabs(matrix_get(U, k, i));
        linha_max = k;
      }
    }

    // Se o maior elemento for zero, a matriz é singular (não tem solução única)
    if (max_val == 0.0) {
      free(P);
      matrix_destroy(U);
      matrix_destroy(L);
      errno = EMATRIXSINGULAR;
      return NULL;
    }

    // Se a linha com o maior elemento não for a linha atual, troca-as
    if (linha_max != i) {
      // Troca as linhas na matriz U
      for (j = 0; j < a->ncols; j++) {
        aux = matrix_get(U, i, j);
        matrix_set(U, i, j, matrix_get(U, linha_max, j));
        matrix_set(U, linha_max, j, aux);
      }
      // Troca as linhas na matriz L (apenas os elementos já calculados antes da coluna i)
      for (j = 0; j < i; j++) {
        aux = matrix_get(L, i, j);
        matrix_set(L, i, j, matrix_get(L, linha_max, j));
        matrix_set(L, linha_max, j, aux);
      }
      // Registra a troca no vetor de permutação P
      int aux_p = P[i];
      P[i] = P[linha_max];
      P[linha_max] = aux_p;
    }
    // -------------------------------------

    // --- PASSO DE ELIMINAÇÃO (Escalonamento) ---
    for (j = i + 1; j < a->nlins; j++) {
      fator = matrix_get(U, j, i) / matrix_get(U, i, i);
      matrix_set(L, j, i, fator); // Guarda o multiplicador na matriz L
      for (k = i; k < a->ncols; k++) {
        matrix_set(U, j, k, matrix_get(U, j, k) - (fator * matrix_get(U, i, k)));
      }
    }
  }

  // 3. Montar a matriz de permutação, que faz parte do resultado: P*A = L*U
  Matrix* permutacao = lu_create_permutation(a->nlins, P);
  free(P);
  if (permutacao == NULL) {
    matrix_destroy(U);
    matrix_destroy(L);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  LU* lu = (LU*)malloc(sizeof(LU));
  if (lu == NULL) {
    matrix_destroy(permutacao);
    matrix_destroy(U);
    matrix_destroy(L);
    errno = EMATRIXALLOCATION;
    return NULL;
  }
  lu->l = L;
  lu->u = U;
  lu->p = permutacao;

  return lu;
}

void lu_destroy(LU* lu) {
  if (lu != NULL) {
    if (lu->l != NULL) {
      matrix_destroy(lu->l);
    }
    if (lu->u != NULL) {
      matrix_destroy(lu->u);
    }
    if (lu->p != NULL) {
      matrix_destroy(lu->p);
    }
    free(lu);
  }
}
