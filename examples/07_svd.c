/**
 * Exemplo 07 - Decomposicao em valores singulares (SVD)
 *
 * Cenario real: compressao de dados de um sensor. Uma matriz A de medidas
 * (6 x 3) e fatorada como A = U * diag(sigma) * V^T pelo algoritmo de
 * Golub-Reinsch. Descartando os valores singulares pequenos, o mesmo
 * conjunto de dados fica comprimido com pouco erro.
 *
 * Este exemplo mostra as tres etapas peas quais a lib expoe a SVD:
 *   svd_householder_bidiag -> svd_qr_bidiag_iter -> svd_pos_proc,
 *   ou o atalho svd_golub_reinsch que encadeia tudo.
 */

#include <stdio.h>
#include <errno.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "svd.h"

/** maior valor absoluto da matriz (medida do erro maximo). */
double erro_max(const Matrix* a, const Matrix* b) {
  double max = 0.0;
  for (unsigned int i = 0; i < a->nlins; i++) {
    for (unsigned int j = 0; j < a->ncols; j++) {
      double e = fabs(matrix_get(a, i, j) - matrix_get(b, i, j));
      if (e > max) {
        max = e;
      }
    }
  }
  return max;
}

int main(void) {
  printf("=== 07 - SVD (Golub-Reinsch) ===\n\n");

  // Leituras de um sensor com 2 sinais dominantes + ruido pequeno
  // (as colunas sao as leituras de cada canal)
  double c1[] = { 3.0, 5.0, -1.0, 8.0, 4.0, 2.0 };
  double c2[] = { 2.0, 4.0, -2.0, 6.0, 3.0, 1.0 };
  double c3[] = { 0.1, 0.2, 0.0, 0.3, 0.1, 0.0 };
  Vector colunas[] = {
    { .dim = 6, .x = c1 },
    { .dim = 6, .x = c2 },
    { .dim = 6, .x = c3 },
  };
  Matrix* a = matrix_create(6, 3, colunas);
  matrix_println(a, "A (6 s leituras x 3 canais)");

  // ---- Atalho: SVD completa em uma chamada ----
  SVDResult* r = svd_golub_reinsch(a);
  if (r == NULL) {
    printf("Falha na SVD (errno=%d)\n", errno);
    return 1;
  }

  matrix_println(r->U, "U (6x3)");
  vector_println(r->sigma, "sigma (valores singulares decrescentes)");
  matrix_println(r->VT, "V^T (3x3)");

  // A = U * diag(sigma) * V^T
  Matrix* sigma = matrix_create_empty(3, 3);
  for (unsigned int i = 0; i < 3; i++) {
    matrix_set(sigma, i, i, r->sigma->x[i]);
  }
  Matrix* us = matrix_mul(r->U, sigma);
  Matrix* usvt = matrix_mul(us, r->VT);
  matrix_println(usvt, "U*diag(sigma)*V^T (espera ~A)");
  printf("Erro maximo da reconstrucao: %.2e\n\n", erro_max(a, usvt));

  // U e V^T sao ortogonais: U^T*U = I e VT*VT^T = I
  Matrix* ut = matrix_transpose(r->U);
  Matrix* utu = matrix_mul(ut, r->U);
  printf("Erro de ortogonalidade de U (max |U^T*U - I|): %.2e\n",
         erro_max(utu, matrix_create_identity(3, 3)));
  matrix_destroy(ut);
  matrix_destroy(utu);
  matrix_destroy(sigma);
  matrix_destroy(us);
  matrix_destroy(usvt);

  // A mesma decomposicao tambem diagnostica o "posto" efetivo: o sigma
  // dois cai ~10x e os demais caem muito mais
  for (unsigned int i = 0; i < r->sigma->dim; i++) {
    printf("sigma[%u] = %.6f\n", i, r->sigma->x[i]);
  }
  svd_destroy(r);

  // ---- Passo a passo, com as etapas individuais ----
  printf("\n--- Etapas individuais ---\n\n");

  HouseholderBidiagResult* hbdr = svd_householder_bidiag(a);
  if (hbdr == NULL) {
    printf("Falha na bidiagonalizacao (errno=%d)\n", errno);
    return 1;
  }
  vector_println(hbdr->d, "diagonal d da bidiagonal B");
  vector_println(hbdr->e, "superdiagonal e da bidiagonal B");
  printf("A = U * B * V^T com B = diag(d) + superdiag(e)\n\n");

  // A iteracao assume a-posse de d, U e V de hbdr; liberar hbdr depois
  QRBidiagResult* bdr = svd_qr_bidiag_iter(hbdr);
  if (bdr == NULL) {
    printf("Falha na iteracao QR (errno=%d)\n", errno);
    svd_householder_bidiag_destroy(hbdr);
    return 1;
  }
  svd_householder_bidiag_destroy(hbdr);

  SVDResult* r2 = svd_pos_proc(bdr);
  if (r2 == NULL) {
    printf("Falha no pos-processamento (errno=%d)\n", errno);
    svd_qr_bidiag_destroy(bdr);
    return 1;
  }
  svd_qr_bidiag_destroy(bdr);

  vector_println(r2->sigma, "sigma (pos-processado)");

  // Rotacoes de Givens: utilzadas durante a iteracao de QR para zerar a
  // superdiagonal
  double c, s;
  svd_givens_rotation(3.0, 4.0, &c, &s);
  printf("Rotacao de Givens p/ (3,4): c=%.4f s=%.4f "
         "(gira (3,4) para (5,0))\n", c, s);
  printf("Sinal (funcao auxiliar): lla_signal(-2.0) = %d, "
         "lla_signal(1.5) = %d\n", lla_signal(-2.0), lla_signal(1.5));

  svd_destroy(r2);
  matrix_destroy(a);

  return 0;
}
