#include "lu.h"

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
    // P(i, p[i]) = 1  ->  em column-major: data[p[i] * n + i]
    P->data[(size_t)p[i] * n + i] = 1.0;
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

  const unsigned int n = a->nlins;
  const size_t total = ((size_t)n * (size_t)n) * sizeof(double);

  Matrix* L = matrix_create_empty(n, n);
  if (L == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  Matrix* U = matrix_create_empty(n, n);
  if (U == NULL) {
    matrix_destroy(L);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  // Vetor de permutação que registra as trocas de linha
  int* P = (int*)malloc(n * sizeof(int));
  if (P == NULL) {
    matrix_destroy(U);
    matrix_destroy(L);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  // Vetor auxiliar com os multiplicadores da coluna corrente
  double* fator = (double*)malloc(n * sizeof(double));
  if (fator == NULL) {
    free(P);
    matrix_destroy(U);
    matrix_destroy(L);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  double* lu_data = U->data;
  double* ll_data = L->data;

  // 1. Copiar a matriz A para U (o escalonamento acontece em U) e preparar
  //    o vetor de permutação e a diagonal unitária de L
  memcpy(lu_data, a->data, total);
  for (unsigned int i = 0; i < n; i++) {
    P[i] = (int)i;
    ll_data[(size_t)i * n + i] = 1.0;
  }

  // 2. Processo de Decomposição LU com Pivoteamento Parcial.
  //    Com column-major, a coluna i de U é o trecho contíguo
  //    U->data + i*n .. + i*n + n-1: a busca de pivô e os multiplicadores
  //    varrem memória contígua, e a atualização do escalonamento corre
  //    por colunas (cada coluna k é atualizada com um laço contíguo).
  for (unsigned int i = 0; i < n; i++) {
    double* col_i = lu_data + (size_t)i * n;

    // --- PASSO DE PIVOTEAMENTO PARCIAL ---
    unsigned int linha_max = i;
    double max_val = fabs(col_i[i]);

    // Procura pelo maior elemento em módulo na coluna i, da linha i para baixo
    for (unsigned int k = i + 1; k < n; k++) {
      if (fabs(col_i[k]) > max_val) {
        max_val = fabs(col_i[k]);
        linha_max = k;
      }
    }

    // Se o maior elemento for zero, a matriz é singular (não tem solução única)
    if (max_val == 0.0) {
      free(fator);
      free(P);
      matrix_destroy(U);
      matrix_destroy(L);
      errno = EMATRIXSINGULAR;
      return NULL;
    }

    // Se a linha com o maior elemento não for a linha atual, troca-as
    if (linha_max != i) {
      // Troca as linhas nas matrizes U e L (as linhas têm passo n num
      // bloque densos column-major, então a troca passa por todas as colunas)
      for (unsigned int j = 0; j < n; j++) {
        const size_t base = (size_t)j * n;
        const double aux_u = lu_data[base + i];
        lu_data[base + i] = lu_data[base + linha_max];
        lu_data[base + linha_max] = aux_u;

        if (j < i) {
          const double aux_l = ll_data[base + i];
          ll_data[base + i] = ll_data[base + linha_max];
          ll_data[base + linha_max] = aux_l;
        }
      }
      // Recarrega a coluna do pivô (o ponteiro pode ter elementos trocados)
      col_i = lu_data + (size_t)i * n;

      // Registra a troca no vetor de permutação P
      int aux_p = P[i];
      P[i] = P[linha_max];
      P[linha_max] = aux_p;
    }
    // -------------------------------------

    // --- PASSO DE ELIMINAÇÃO (Escalonamento) ---
    // Multiplicadores da coluna i: L[j][i] = U[j][i] / U[i][i]
    const double piv_val = col_i[i];
    for (unsigned int j = i + 1; j < n; j++) {
      const double f = col_i[j] / piv_val;
      fator[j] = f;
      ll_data[(size_t)i * n + j] = f;
      col_i[j] = 0.0; // zera a coluna de U abaixo do pivô de forma definitiva
    }

    // Atualiza o restante de U por colunas: para cada coluna k de i+1 até
    // n-1, U[j][k] -= fator[j] * U[i][k], com a coluna acessada contiguamente.
    for (unsigned int k = i + 1; k < n; k++) {
      double* col_k = lu_data + (size_t)k * n;
      const double uik = col_k[i];
      if (uik == 0.0) {
        continue;
      }
      for (unsigned int j = i + 1; j < n; j++) {
        col_k[j] -= fator[j] * uik;
      }
    }
  }

  // 3. Montar a matriz de permutação, que faz parte do resultado: P*A = L*U.
  //    P só é liberado no erro ou transferido para o resultado (lu->piv);
  //    nunca após um free.
  Matrix* permutacao = lu_create_permutation(n, P);
  if (permutacao == NULL) {
    free(P);
    free(fator);
    matrix_destroy(U);
    matrix_destroy(L);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  LU* lu = (LU*)malloc(sizeof(LU));
  if (lu == NULL) {
    free(P);
    free(fator);
    matrix_destroy(permutacao);
    matrix_destroy(U);
    matrix_destroy(L);
    errno = EMATRIXALLOCATION;
    return NULL;
  }
  lu->l = L;
  lu->u = U;
  lu->p = permutacao;
  lu->piv = P; // o vetor de trocas de linha faz parte do resultado agora

  free(fator);

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
    if (lu->piv != NULL) {
      free(lu->piv);
    }
    free(lu);
  }
}
