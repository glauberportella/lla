#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include "vector.h"
#include "matrix.h"
#include "lu.h"

void test_vector_create_empty() {
  printf("--- TEST: vector_create_empty ---\n");

  Vector* a = vector_create_empty(3);
  assert(a != NULL);
  assert(a->dim == 3);

  for (unsigned int i = 0; i < a->dim; i++) {
    assert(a->x[i] == 0.0);
  }

  vector_destroy(a);

  // Dimensao invalida
  assert(vector_create_empty(0) == NULL);
  assert(errno == EINVALIDDIM);
}

void test_vector_create() {
  printf("--- TEST: vector_create ---\n");

  double coords[] = { 1.0, -1.0 };
  Vector* a = vector_create(2, coords);
  assert(a != NULL);
  assert(a->dim == 2);
  assert(a->x[0] == 1.0);
  assert(a->x[1] == -1.0);

  vector_destroy(a);

  // Dimensao invalida
  double coords_invalidos[] = { 1.0 };
  assert(vector_create(0, coords_invalidos) == NULL);
  assert(errno == EINVALIDDIM);
}

void test_vector_add() {
  printf("--- TEST: vector_add ---\n");

  double coords_a[] = { 1.0, -1.0 };
  Vector* a = vector_create(2, coords_a);
  assert(a != NULL);

  double coords_b[] = { -1.0, 1.0 };
  Vector* b = vector_create(2, coords_b);
  assert(b != NULL);

  Vector* c = vector_add(a, b);
  assert(c != NULL);
  assert(c->dim == 2);
  assert(c->x[0] == 0.0);
  assert(c->x[1] == 0.0);

  // Vetor nulo
  assert(vector_add(a, NULL) == NULL);
  assert(errno == EVECTORNULL);

  // Dimensao invalida
  double coords_d[] = { 1.0, 2.0, 3.0 };
  Vector* d = vector_create(3, coords_d);
  assert(d != NULL);

  assert(vector_add(a, d) == NULL);
  assert(errno == EINVALIDDIM);

  vector_destroy(d);
  vector_destroy(c);
  vector_destroy(b);
  vector_destroy(a);
}

void test_vector_sub() {
  printf("--- TEST: vector_sub ---\n");

  double coords_a[] = { 1.0, -1.0 };
  Vector* a = vector_create(2, coords_a);
  assert(a != NULL);

  double coords_b[] = { -1.0, 1.0 };
  Vector* b = vector_create(2, coords_b);
  assert(b != NULL);

  Vector* c = vector_sub(a, b);
  assert(c != NULL);
  assert(c->dim == 2);
  assert(c->x[0] == 2.0);
  assert(c->x[1] == -2.0);

  // Vetor nulo
  assert(vector_sub(a, NULL) == NULL);
  assert(errno == EVECTORNULL);

  // Dimensao invalida
  double coords_d[] = { 1.0, 2.0, 3.0 };
  Vector* d = vector_create(3, coords_d);
  assert(d != NULL);

  assert(vector_sub(a, d) == NULL);
  assert(errno == EINVALIDDIM);

  vector_destroy(d);
  vector_destroy(c);
  vector_destroy(b);
  vector_destroy(a);
}

void test_vector_mul() {
  printf("--- TEST: vector_mul ---\n");

  double coords[] = { -1.0, 2.0 };
  Vector* a = vector_create(2, coords);
  assert(a != NULL);
  vector_println(a, "A");

  Vector* b = vector_mul(-2.0, a);
  assert(b != NULL);
  assert(b->dim == 2);
  assert(b->x[0] == 2.0);
  assert(b->x[1] == -4.0);
  vector_println(b, "B = -2A");

  // Escalar zero
  Vector* c = vector_mul(0.0, a);
  assert(c != NULL);
  assert(c->x[0] == 0.0);
  assert(c->x[1] == 0.0);

  // Vetor nulo
  assert(vector_mul(2.0, NULL) == NULL);
  assert(errno == EVECTORNULL);

  vector_destroy(c);
  vector_destroy(b);
  vector_destroy(a);
}

void test_vector_inner_product() {
  printf("--- TEST: vector_inner_product ---\n");

  double coords_a[] = { 1.0, 2.0, 3.0 };
  Vector* a = vector_create(3, coords_a);
  assert(a != NULL);

  double coords_b[] = { 4.0, 5.0, 6.0 };
  Vector* b = vector_create(3, coords_b);
  assert(b != NULL);

  // 1*4 + 2*5 + 3*6 = 32
  double p = vector_inner_product(a, b);
  assert(!isnan(p));
  assert(fabs(p - 32.0) < 1e-9);
  printf("A . B = %.2f\n", p);

  // Vetores ortogonais
  double coords_c[] = { 1.0, 0.0, 0.0 };
  Vector* c = vector_create(3, coords_c);
  assert(c != NULL);

  double coords_d[] = { 0.0, 1.0, 0.0 };
  Vector* d = vector_create(3, coords_d);
  assert(d != NULL);

  assert(fabs(vector_inner_product(c, d)) < 1e-9);

  // Dimensao invalida
  double coords_e[] = { 1.0, 2.0 };
  Vector* e = vector_create(2, coords_e);
  assert(e != NULL);

  assert(isnan(vector_inner_product(a, e)));
  assert(errno == EINVALIDDIM);

  // Vetor nulo
  assert(isnan(vector_inner_product(a, NULL)));
  assert(errno == EVECTORNULL);

  vector_destroy(e);
  vector_destroy(d);
  vector_destroy(c);
  vector_destroy(b);
  vector_destroy(a);
}

void test_vector_norm_p() {
  printf("--- TEST: vector_norm_p ---\n");

  double coords[] = { -3.0, 4.0 };
  Vector* a = vector_create(2, coords);
  assert(a != NULL);

  // p = 1: 3 + 4 = 7
  double n1 = vector_norm_p(a, 1);
  assert(!isnan(n1));
  assert(fabs(n1 - 7.0) < 1e-9);
  printf("||A||_1 = %.2f\n", n1);

  // p = 2: sqrt(9 + 16) = 5
  double n2 = vector_norm_p(a, 2);
  assert(!isnan(n2));
  assert(fabs(n2 - 5.0) < 1e-9);
  printf("||A||_2 = %.2f\n", n2);

  // p = 3: (27 + 64)^(1/3)
  double n3 = vector_norm_p(a, 3);
  assert(!isnan(n3));
  assert(fabs(n3 - 4.497941445275415) < 1e-9);
  printf("||A||_3 = %.2f\n", n3);

  // Vetor nulo
  assert(isnan(vector_norm_p(NULL, 2)));
  assert(errno == EVECTORNULL);

  // Ordem invalida
  assert(isnan(vector_norm_p(a, 0)));
  assert(errno == EINVALIDP);

  vector_destroy(a);
}

void test_vector_norm_inf() {
  printf("--- TEST: vector_norm_inf ---\n");

  double coords[] = { -4.0, 3.0 };
  Vector* a = vector_create(2, coords);
  assert(a != NULL);

  // max(|x_i|) = 4
  double n = vector_norm_inf(a);
  assert(!isnan(n));
  assert(fabs(n - 4.0) < 1e-9);
  printf("||A||_inf = %.2f\n", n);

  // Todas as coordenadas negativas
  double coords_b[] = { -4.0, -3.0 };
  Vector* b = vector_create(2, coords_b);
  assert(b != NULL);

  double n_b = vector_norm_inf(b);
  assert(!isnan(n_b));
  assert(fabs(n_b - 4.0) < 1e-9);

  // Vetor zerado
  double coords_c[] = { 0.0, 0.0 };
  Vector* c = vector_create(2, coords_c);
  assert(c != NULL);

  assert(fabs(vector_norm_inf(c)) < 1e-9);

  // Vetor nulo passado como argumento
  assert(isnan(vector_norm_inf(NULL)));
  assert(errno == EVECTORNULL);

  vector_destroy(c);
  vector_destroy(b);
  vector_destroy(a);
}

void test_vector_norm() {
  printf("--- TEST: vector_norm ---\n");

  double coords[] = { 3.0, 4.0 };
  Vector* a = vector_create(2, coords);
  assert(a != NULL);

  // sqrt(9 + 16) = 5
  double n = vector_norm(a);
  assert(!isnan(n));
  assert(fabs(n - 5.0) < 1e-9);
  printf("||A|| = %.2f\n", n);

  // Equivale a vector_norm_p com p = 2
  assert(fabs(vector_norm(a) - vector_norm_p(a, 2)) < 1e-9);

  // Vetor nulo
  assert(isnan(vector_norm(NULL)));
  assert(errno == EVECTORNULL);

  vector_destroy(a);
}

void test_vector_print() {
  printf("--- TEST: vector_print ---\n");

  double coords_a[] = { 1.0, -1.0 };
  Vector* a = vector_create(2, coords_a);
  assert(a != NULL);
  vector_println(a, "A");

  double coords_b[] = { -1.0, 1.0 };
  Vector* b = vector_create(2, coords_b);
  assert(b != NULL);
  vector_println(b, "B");

  Vector* c = vector_add(a, b);
  assert(c != NULL);
  vector_println(c, "C = A + B");

  Vector* d = vector_sub(a, b);
  assert(d != NULL);
  vector_println(d, "D = A - B");

  // Sem rotulo
  vector_println(a, NULL);

  vector_destroy(d);
  vector_destroy(c);
  vector_destroy(b);
  vector_destroy(a);
}

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

  // Matriz que exige pivotamento: det = -1
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
  test_vector_create_empty();
  test_vector_create();
  test_vector_add();
  test_vector_sub();
  test_vector_mul();
  test_vector_inner_product();
  test_vector_norm_p();
  test_vector_norm_inf();
  test_vector_norm();
  test_vector_print();
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
