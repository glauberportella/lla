#include "vector.h"
#include "matrix.h"
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

Vector* vector_product(unsigned int dim, const Vector vectors[]) {
  if (dim < 2) {
    errno = EINVALIDDIM;
    return NULL;
  }

  // Vamos montar uma matriz com os vetores de entrada sendo que
  //     |i   j   k  ... |
  // M = |a1  a2  a3 ... |
  //     |b1  b2  b3 ... |
  
  /* OPCAO 1
  Matrix* m = matrix_create_empty(dim - 1, dim);
  if (m == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }
  // temos que preencher a matriz com os vetores de entrada como vetores linha
  for (i = 0; i < dim - 1; i++) {
    Vector v = vectors[i];
    for (j = 0; j < dim; j++) {
      matrix_set(m, j, i, v.x[i]);
    }
  }
  */
  // OPCAO 2
  Matrix* temp = matrix_create(dim, dim - 1, (Vector*)vectors);
  if (temp == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }
  Matrix* aux = matrix_transpose(temp);
  if (aux == NULL) {
    matrix_destroy(temp);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  // debug
  double det = 0.00;
  int sinal_cof = 1;
  double* coords = (double*)calloc(dim, sizeof(double));
  if (coords == NULL) {
    matrix_destroy(aux);
    matrix_destroy(temp);
    errno = EVECTORALLOCATION;
    return NULL;
  }

  for (unsigned int i = 0; i < dim; i++) {
    // monta submatriz eliminando a coluna i
    Matrix* submatrix = matrix_remove_column_vector_at(aux, i);
    if (submatrix == NULL) {
      free(coords);
      matrix_destroy(aux);
      matrix_destroy(temp);
      errno = EMATRIXALLOCATION;
      return NULL;
    }
    matrix_println(submatrix, "m");
    // calcular determinante da submatriz
    det = matrix_determinant(submatrix);
    sinal_cof = ((i + 1) + 1 % 2 == 0) ? 1 : -1;
    coords[i] = sinal_cof * det;
    matrix_destroy(submatrix);
  }

  matrix_destroy(aux);
  matrix_destroy(temp);

  Vector* result = vector_create(dim, coords);
  if (result == NULL) {
    errno = EVECTORALLOCATION;
    return NULL;
  }

  vector_println(result, "a x b");
  return result;
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

