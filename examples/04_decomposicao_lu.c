/**
 * Exemplo 04 - Decomposicao LU
 *
 * Cenario real: analise de circuito eletrico. As equacoes de malha de um
 * circuito com 3 nos resulta em um sistema linear cuja matriz e
 * decomposta em LU com pivoteamento parcial: P*A = L*U.
 *
 * Fatores de lu_create_decomposition:
 *  - l: triangular inferior com diagonal unitaria
 *  - u: triangular superior
 *  - p: matriz de permutacao, tal que p * A = L * U
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "lu.h"
#include "lsystem.h"

double tolerancia = 1e-10;

/** Resolve L * y = Pb (triangular inferior, diagonal unitaria). */
void substituicao_forward(const Matrix* l, const Vector* pb, Vector* y) {
  for (unsigned int i = 0; i < y->dim; i++) {
    double soma = 0.0;
    for (unsigned int j = 0; j < i; j++) {
      soma += matrix_get(l, i, j) * y->x[j];
    }
    y->x[i] = pb->x[i] - soma; // l(i,i) = 1
  }
}

/** Resolve U * x = y (triangular superior). */
void substituicao_back(const Matrix* u, const Vector* y, Vector* x) {
  for (int i = (int)u->ncols - 1; i >= 0; i--) {
    double soma = 0.0;
    for (unsigned int j = (unsigned int)i + 1; j < u->ncols; j++) {
      soma += matrix_get(u, (unsigned int)i, j) * x->x[j];
    }
    x->x[(unsigned int)i] = (y->x[(unsigned int)i] - soma) / matrix_get(u, (unsigned int)i, (unsigned int)i);
  }
}

int comparar_vetores(const Vector* a, const Vector* b, double tol) {
  if (a == NULL || b == NULL || a->dim != b->dim) {
    return 0;
  }
  for (unsigned int i = 0; i < a->dim; i++) {
    if (fabs(a->x[i] - b->x[i]) > tol) {
      return 0;
    }
  }
  return 1;
}

int main(void) {
  printf("=== 04 - Decomposicao LU ===\n\n");

  // Sistema de malhas (condutancias em mS, correntes em mA):
  //  4*I1 - 1*I2 + 0*I3 = 30
  // -1*I1 + 5*I2 - 1*I3 = 0
  //  0*I1 - 1*I2 + 6*I3 = 20
  double c1[] = { 4.0, -1.0, 0.0 };
  double c2[] = { -1.0, 5.0, -1.0 };
  double c3[] = { 0.0, -1.0, 6.0 };
  Vector colunas[] = {
    { .dim = 3, .x = c1 },
    { .dim = 3, .x = c2 },
    { .dim = 3, .x = c3 },
  };
  Matrix* a = matrix_create(3, 3, colunas);

  double b_coords[] = { 30.0, 0.0, 20.0 };
  Vector* b = vector_create(3, b_coords);

  LU* lu = lu_create_decomposition(a);
  if (lu == NULL) {
    printf("Falha na decomposicao LU (errno=%d)\n", errno);
    return 1;
  }

  matrix_println(lu->p, "P");
  matrix_println(lu->l, "L");
  matrix_println(lu->u, "U");

  // Verificacao: P*A = L*U
  Matrix* pa = matrix_mul(lu->p, a);
  Matrix* lu_prod = matrix_mul(lu->l, lu->u);
  matrix_println(lu_prod, "L*U (espera P*A)");

  // Tambem: A = P^t * L * U, pois P^t = P^-1
  Matrix* pt = matrix_transpose(lu->p);
  Matrix* pt_l = matrix_mul(pt, lu->l);
  Matrix* pt_lu = matrix_mul(pt_l, lu->u);
  printf("A - P^t*L*U fica abaixo do limite?(ver abaixo)\n");
  matrix_println(pt_lu, "P^t*L*U");
  printf("\n");

  // Resolucao do sistema via substituicoes sobre L e U:
  // 1) P*b
  Vector pb_struct = { .dim = 3, .x = (double*)malloc(3 * sizeof(double)) };
  for (unsigned int i = 0; i < 3; i++) { // P*b
    double soma = 0.0;
    for (unsigned int j = 0; j < 3; j++) {
      soma += matrix_get(lu->p, i, j) * b->x[j];
    }
    pb_struct.x[i] = soma;
  }

  Vector* y = vector_create_empty(3);
  substituicao_forward(lu->l, &pb_struct, y);
  vector_println(y, "L*y = P*b -> y");

  Vector* x = vector_create_empty(3);
  substituicao_back(lu->u, y, x);
  vector_println(x, "U*x = y -> x (correntes I em mA)");

  // Cross-check com lsystem_solve
  Vector* x_check = lsystem_solve(a, b);
  printf("x coincide com lsystem_solve? %s\n",
         comparar_vetores(x, x_check, tolerancia) ? "sim" : "nao");

  vector_destroy(b);
  vector_destroy(y);
  vector_destroy(x);
  vector_destroy(x_check);
  free(pb_struct.x);
  matrix_destroy(a);
  matrix_destroy(pa);
  matrix_destroy(lu_prod);
  matrix_destroy(pt);
  matrix_destroy(pt_l);
  matrix_destroy(pt_lu);
  lu_destroy(lu);

  return 0;
}
