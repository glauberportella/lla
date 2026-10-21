#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include "vector.h"
#include "matrix.h"
#include "lu.h"

void test_matrix_create() {
  printf("--- TEST: matrix_create ---\n");

  double v1[] = { 1.0, -1.0 };
  double v2[] = { -1.0, 1.0 };
  Vector vetor_colunas[] = {
    { .dim = 2, .x = v1 },
    { .dim = 2, .x = v2 },
  };
  Matrix* a = matrix_create(2, 2, vetor_colunas);

  assert(a != NULL);
  assert(a->nlins == 2);
  assert(a->ncols == 2);
  assert(a->columns[0].dim == 2);
  assert(a->columns[1].dim == 2);
  assert(a->columns[0].x[0] == 1.0);
  assert(a->columns[0].x[1] == -1.0);
  assert(a->columns[1].x[0] == -1.0);
  assert(a->columns[1].x[1] == 1.0);

  matrix_destroy(a);
}

void test_matrix_add() {
  printf("--- TEST: matrix_add ---\n");

  double v1[] = { 1.0, -1.0 };
  double v2[] = { -1.0, 1.0 };
  Vector vetor_colunas_a[] = {
    { .dim = 2, .x = v1 },
    { .dim = 2, .x = v2 },
  };
  Matrix* a = matrix_create(2, 2, vetor_colunas_a);
  assert(a != NULL);

  double v3[] = { -1.0, 0.0 };
  double v4[] = { 0.0, -1.0 };
  Vector vetor_colunas_b[] = {
    { .dim = 2, .x = v3 },
    { .dim = 2, .x = v4 },
  };
  Matrix* b = matrix_create(2, 2, vetor_colunas_b);
  assert(b != NULL);

  Matrix* c = matrix_add(a, b);
  assert(c != NULL);

  assert(c->columns[0].x[0] == 0.0);
  assert(c->columns[0].x[1] == -1.0);
  assert(c->columns[1].x[0] == -1.0);
  assert(c->columns[1].x[1] == 0.0);

  matrix_destroy(c);
  matrix_destroy(b);
  matrix_destroy(a);
}

void test_matrix_sub() {
  printf("--- TEST: matrix_sub ---\n");

  double v1[] = { 1.0, -1.0 };
  double v2[] = { -1.0, 1.0 };
  Vector vetor_colunas_a[] = {
    { .dim = 2, .x = v1 },
    { .dim = 2, .x = v2 },
  };
  Matrix* a = matrix_create(2, 2, vetor_colunas_a);
  assert(a != NULL);

  double v3[] = { -1.0, 0.0 };
  double v4[] = { 0.0, -1.0 };
  Vector vetor_colunas_b[] = {
    { .dim = 2, .x = v3 },
    { .dim = 2, .x = v4 },
  };
  Matrix* b = matrix_create(2, 2, vetor_colunas_b);
  assert(b != NULL);

  Matrix* c = matrix_sub(a, b);
  assert(c != NULL);

  assert(c->columns[0].x[0] == 2.0);
  assert(c->columns[0].x[1] == -1.0);
  assert(c->columns[1].x[0] == -1.0);
  assert(c->columns[1].x[1] == 2.0);

  matrix_destroy(c);
  matrix_destroy(b);
  matrix_destroy(a);
}

void test_matrix_mul_scalar() {
  printf("--- TEST: matrix_mul_scalar ---\n");

  double v1[] = { -1.0, -1.0 };
  double v2[] = { -1.0, -1.0 };
  Vector vetor_colunas[] = {
    { .dim = 2, .x = v1 },
    { .dim = 2, .x = v2 },
  };
  
  Matrix* a = matrix_create(2, 2, vetor_colunas);
  assert(a != NULL);
  
  matrix_println(a, "A");

  Matrix* b = matrix_mul_scalar(-2.0, a);
  assert(b != NULL);
  
  assert(b->columns[0].x[0] == 2.0);
  assert(b->columns[0].x[1] == 2.0);
  assert(b->columns[1].x[0] == 2.0);
  assert(b->columns[1].x[1] == 2.0);

  matrix_println(b, "B = -2A");

  matrix_destroy(b);
  matrix_destroy(a);
}

void test_matrix_mul() {
  printf("--- TEST: matrix_mul ---\n");

  double v1[] = { -1.0, -1.0, 1.0 };
  double v2[] = { -1.0, -1.0, 1.0 };
  Vector vetor_colunas_a[] = {
    { .dim = 3, .x = v1 },
    { .dim = 3, .x = v2 },
  };
  
  Matrix* a = matrix_create(3, 2, vetor_colunas_a);
  assert(a != NULL);
  matrix_println(a, "A");

  double v3[] = { 1.0, 0.0 };
  double v4[] = { 0.0, -1.0 };
  double v5[] = { -1.0, -1.0 };
  Vector vetor_colunas_b[] = {
    { .dim = 2, .x = v3 },
    { .dim = 2, .x = v4 },
    { .dim = 2, .x = v5 },
  };

  Matrix* b = matrix_create(2, 3, vetor_colunas_b);
  assert(b != NULL);
  matrix_println(b, "B");

  Matrix* ab = matrix_mul(a, b);
  assert(ab != NULL);
  matrix_println(ab, "AxB");

  Matrix* ba = matrix_mul(b, a);
  assert(ba != NULL);
  matrix_println(ba, "BxA");

  matrix_destroy(ba);
  matrix_destroy(ab);
  matrix_destroy(b);
  matrix_destroy(a);
}

void test_matrix_print() {
  printf("--- TEST: matrix_print ---\n");

  double v1[] = { 1.0, -1.0 };
  double v2[] = { -1.0, 1.0 };
  Vector vetor_colunas[] = {
    { .dim = 2, .x = v1 },
    { .dim = 2, .x = v2 },
  };
  
  Matrix* a = matrix_create(2, 2, vetor_colunas);
  assert(a != NULL);
  matrix_println(a, "A");

  double v3[] = { -1.0, 0.0 };
  double v4[] = { 0.0, -1.0 };
  Vector vetor_colunas_b[] = {
    { .dim = 2, .x = v3 },
    { .dim = 2, .x = v4 },
  };
  Matrix* b = matrix_create(2, 2, vetor_colunas_b);
  assert(b != NULL);
  matrix_println(b, "B");
  
  Matrix* c = matrix_add(a, b);
  assert(c != NULL);
  matrix_println(c, "C = A + B");

  Matrix* d = matrix_sub(a, b);
  assert(d != NULL);
  matrix_println(d, "D = A - B");

  matrix_destroy(d);
  matrix_destroy(c);
  matrix_destroy(b);
  matrix_destroy(a);
}

void test_matrix_lu_decomposition() {
  printf("--- TEST: matrix_lu_decomposition ---\n");

  double v1[] = {2.0, 4.0, 8.0};
  double v2[] = {1.0, 3.0, 7.0};
  double v3[] = {1.0, 1.0, 3.0};
  Vector vetor_colunas[] = {
    { .dim = 3, .x = v1 },
    { .dim = 3, .x = v2 },
    { .dim = 3, .x = v3 },
  };
  
  Matrix* a = matrix_create(3, 3, vetor_colunas);
  assert(a != NULL);
  matrix_println(a, "A");

  LU* lu = lu_create_decomposition(a);
  if (lu == NULL) {
    printf("ERRO: %d", errno);
  }
  assert(lu != NULL);

  matrix_println(lu->l, "L");
  matrix_println(lu->u, "U");
  matrix_println(lu->p, "P");

  // P*A = L*U
  Matrix* pa = matrix_mul(lu->p, a);
  Matrix* lu_prod = matrix_mul(lu->l, lu->u);
  assert(pa != NULL);
  assert(lu_prod != NULL);

  for (unsigned int i = 0; i < a->nlins; i++) {
    for (unsigned int j = 0; j < a->ncols; j++) {
      assert(fabs(matrix_get(pa, i, j) - matrix_get(lu_prod, i, j)) < 1e-9);
    }
  }

  matrix_destroy(lu_prod);
  matrix_destroy(pa);
  lu_destroy(lu);
  matrix_destroy(a);
}

void test_matrix_transpose() {
  printf("--- TEST: matrix_transpose ---\n");

  double v1[] = {2.0, 4.0, 8.0};
  double v2[] = {1.0, 3.0, 7.0};
  double v3[] = {1.0, 1.0, 3.0};
  Vector vetor_colunas[] = {
    { .dim = 3, .x = v1 },
    { .dim = 3, .x = v2 },
    { .dim = 3, .x = v3 },
  };
  
  Matrix* a = matrix_create(3, 3, vetor_colunas);
  assert(a != NULL);
  matrix_println(a, "A");

  Matrix* at = matrix_transpose(a);
  assert(at != NULL);
  matrix_println(at, "At");

  assert(at->nlins == a->ncols);
  assert(at->ncols == a->nlins);

  for (unsigned int i = 0; i < a->nlins; i++) {
    for (unsigned int j = 0; j < a->ncols; j++) {
      assert(fabs(matrix_get(at, i, j) - matrix_get(a, j, i)) < 1e-9);
    }
  }

  matrix_destroy(at);
  matrix_destroy(a);
}

void test_matrix_inverse() {
  printf("--- TEST: matrix_inverse ---\n");

  double v1[] = {2.0, 4.0, 8.0};
  double v2[] = {1.0, 3.0, 7.0};
  double v3[] = {1.0, 1.0, 3.0};
  Vector vetor_colunas[] = {
    { .dim = 3, .x = v1 },
    { .dim = 3, .x = v2 },
    { .dim = 3, .x = v3 },
  };
  
  Matrix* a = matrix_create(3, 3, vetor_colunas);
  assert(a != NULL);
  matrix_println(a, "A");

  Matrix* inverse = matrix_inverse(a);
  assert(inverse != NULL);
  matrix_println(inverse, "A^-1");

  // A * A^-1 = I
  Matrix* product = matrix_mul(a, inverse);
  assert(product != NULL);

  for (unsigned int i = 0; i < a->nlins; i++) {
    for (unsigned int j = 0; j < a->ncols; j++) {
      double expected = (i == j) ? 1.0 : 0.0;
      assert(fabs(matrix_get(product, i, j) - expected) < 1e-9);
    }
  }

  matrix_destroy(product);
  matrix_destroy(inverse);
  matrix_destroy(a);
}

void test_matrix_determinant() {
  printf("--- TEST: matrix_determinant ---\n");

  double v1[] = {-3.0, 1.0, -2.0};
  double v2[] = {2.0, 5.0, 1.5};
  double v3[] = {3.0, 2.0, 2.0};
  Vector vetor_colunas[] = {
    { .dim = 3, .x = v1 },
    { .dim = 3, .x = v2 },
    { .dim = 3, .x = v3 },
  };
  
  Matrix* a = matrix_create(3, 3, vetor_colunas);
  assert(a != NULL);
  matrix_println(a, "A");

  // det(A) = 1.5
  double det = matrix_determinant(a);
  assert(!isnan(det));
  assert(fabs(det - 1.5) < 1e-9);
  printf("Det(A) = %.2f\n", det);

  matrix_destroy(a);

  // Matriz que exige pivoteamento: det = -1
  double p1[] = {0.0, 1.0};
  double p2[] = {1.0, 0.0};
  Vector colunas_piv[] = {
    { .dim = 2, .x = p1 },
    { .dim = 2, .x = p2 },
  };

  Matrix* b = matrix_create(2, 2, colunas_piv);
  assert(b != NULL);

  double det_b = matrix_determinant(b);
  assert(!isnan(det_b));
  assert(fabs(det_b + 1.0) < 1e-9);
  printf("Det(B) = %.2f\n", det_b);

  matrix_destroy(b);

  // Matriz singular: det = 0
  double s1[] = {1.0, 2.0};
  double s2[] = {2.0, 4.0};
  Vector colunas_sing[] = {
    { .dim = 2, .x = s1 },
    { .dim = 2, .x = s2 },
  };

  Matrix* c = matrix_create(2, 2, colunas_sing);
  assert(c != NULL);

  double det_c = matrix_determinant(c);
  assert(fabs(det_c) < 1e-9);
  printf("Det(C) = %.2f\n", det_c);

  matrix_destroy(c);

  // Matriz nula e matriz nao quadrada
  assert(isnan(matrix_determinant(NULL)));
  assert(errno == EMATRIXNULL);

  double n1[] = {1.0, 2.0};
  Vector colunas_ret[] = {
    { .dim = 2, .x = n1 },
  };

  Matrix* d = matrix_create(2, 1, colunas_ret);
  assert(d != NULL);

  assert(isnan(matrix_determinant(d)));
  assert(errno == EMATRIXINVALIDDIM);

  matrix_destroy(d);
}

int main() {
  test_matrix_create();
  test_matrix_add();
  test_matrix_sub();
  test_matrix_mul();
  test_matrix_mul_scalar();
  test_matrix_print();
  test_matrix_lu_decomposition();
  test_matrix_transpose();
  test_matrix_inverse();
  test_matrix_determinant();
  return 0;
}
