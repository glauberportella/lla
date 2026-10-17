#include <stdio.h>

#include <assert.h>
#include "vector.h"
#include "matrix.h"

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


int main() {
  test_matrix_create();
  test_matrix_add();
  test_matrix_sub();
  test_matrix_mul();
  test_matrix_mul_scalar();
  test_matrix_print();
  return 0;
}
