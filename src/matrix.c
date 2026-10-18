#include "matrix.h"

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

Matrix* matrix_create_empty(unsigned int nlins, unsigned int ncols) {
  if (nlins == 0 || ncols == 0) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  Matrix* m = (Matrix*)malloc(sizeof(Matrix));
  if (m == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  m->nlins = nlins;
  m->ncols = ncols;
  m->columns = (Vector*)calloc(ncols, sizeof(Vector));

  if (m->columns == NULL) {
    free(m);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  for (unsigned int j = 0; j < ncols; j++) {
    m->columns[j].dim = nlins;
    m->columns[j].x = (double*)calloc(nlins, sizeof(double));

    if (m->columns[j].x == NULL) {
      for (unsigned int k = 0; k < j; k++) {
        free(m->columns[k].x);
      }
      free(m->columns);
      free(m);
      errno = EMATRIXALLOCATION;
      return NULL;
    }
  }

  return m;
}

Matrix* matrix_create(unsigned int nlins, unsigned int ncols, Vector columns[]) {
  if (ncols == 0 || columns == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  Matrix* m = (Matrix*)malloc(sizeof(Matrix));
  if (m == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  m->columns = (Vector*)calloc(ncols, sizeof(Vector));
  if (m->columns == NULL) {
    free(m);
    errno = EMATRIXALLOCATION;
    return NULL;
  }
  m->nlins = nlins;
  m->ncols = ncols;

  for (unsigned int i = 0; i < ncols; i++) {
    if (columns[i].dim != m->nlins || columns[i].x == NULL) {
      for (unsigned int j = 0; j < i; j++) {
        free(m->columns[i].x);
      }
      free(m->columns);
      free(m);
      errno = EMATRIXALLOCATION;
      return NULL;
    }

    m->columns[i].dim = columns[i].dim;
    m->columns[i].x = (double*)malloc(m->columns[i].dim * sizeof(double));
    if (m->columns[i].x == NULL) {
      free(m->columns);
      free(m);
      errno = EMATRIXALLOCATION;
      return NULL;
    }

    memcpy(m->columns[i].x, columns[i].x, columns[i].dim * sizeof(double));
  }

  return m;
}

void matrix_destroy(Matrix* m) {
  if (m != NULL) {
    if (m->columns != NULL) {
      for (unsigned int i = 0; i < m->ncols; i++) {
        if (m->columns[i].x != NULL) {
          free(m->columns[i].x);
        }
      }
      free(m->columns);
    }
    free(m);
  }
}

int matrix_validate_index(const Matrix* m, const unsigned int i, const unsigned int j) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return -2;
  }

  if (i > m->nlins) {
    errno = EMATRIXINVALIDLINEINDEX;
    return -1;
  }

  if (j > m->ncols) {
    errno = EMATRIXINVALIDCOLINDEX;
    return 1;
  }

  return 0;
}

Matrix* matrix_set(Matrix* m, const unsigned int i, const unsigned int j, const double aij) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  const int validated = matrix_validate_index(m, i, j);
  if (validated != 0) {
    return NULL;
  }

  Vector* column_vector = &m->columns[j];
  if (column_vector != NULL) {
    column_vector->x[i] = aij;
  }

  return m;
}

double matrix_get(const Matrix* m, const unsigned int i, const unsigned int j) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return NAN;
  }

  const int validated = matrix_validate_index(m, i, j);
  if (validated != 0) {
    return NAN;
  }

  Vector* column_vector = &m->columns[j];
  if (column_vector == NULL) {
    return NAN;
  }

  return column_vector->x[i];
}

Matrix* matrix_add(const Matrix* a, const Matrix* b) {
  if (a == NULL || b == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  if (a->nlins != b->nlins || a->ncols != b->ncols) {
    errno = EMATRIXINVALIDADDDIM;
    return NULL;
  }

  Matrix* c = matrix_create(a->nlins, a->ncols, a->columns);
  if (c == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  for (unsigned int i = 0; i < a->nlins; i++) {
    for (unsigned int j = 0; j < a->ncols; j++) {
      double aij = matrix_get(a, i, j);
      double bij = matrix_get(b, i, j);
      double cij = aij + bij;
      matrix_set(c, i, j, cij);
    }
  }

  return c;
}

Matrix* matrix_sub(const Matrix* a, const Matrix* b) {
  if (a == NULL || b == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  if (a->nlins != b->nlins || a->ncols != b->ncols) {
    errno = EMATRIXINVALIDADDDIM;
    return NULL;
  }

  Matrix* c = matrix_create(a->nlins, a->ncols, a->columns);
  if (c == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  for (unsigned int i = 0; i < a->nlins; i++) {
    for (unsigned int j = 0; j < a->ncols; j++) {
      double aij = matrix_get(a, i, j);
      double bij = matrix_get(b, i, j);
      matrix_set(c, i, j, aij - bij);
    }
  }

  return c;
}

Matrix* matrix_mul(const Matrix* a, const Matrix* b) {
  if (a == NULL || b == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  if (a->ncols != b->nlins) {
    errno = EMATRIXINVALIDMUL;
    return NULL;
  }

  Matrix* c = matrix_create_empty(a->nlins, b->ncols);
  if (c == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  for (unsigned int i = 0; i < a->nlins; i++) {
    for (unsigned int j = 0; j < b->ncols; j++) {
      double cij = 0.0;
      for (unsigned int k = 0; k < a->ncols; k++) {
        double aik = matrix_get(a, i, k);
        double bkj = matrix_get(b, k, j);
        cij += aik * bkj;
      }
      matrix_set(c, i, j, cij);
    }
  }

  return c;
}

Matrix* matrix_mul_scalar(const double scalar, const Matrix* m) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  Matrix* c = matrix_create(m->nlins, m->ncols, m->columns);
  if (c == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  for (unsigned int i = 0; i < m->nlins; i++) {
    for (unsigned int j = 0; j < m->ncols; j++) {
      double cij = scalar * matrix_get(m, i, j);
      matrix_set(c, i, j, cij);
    }
  }

  return c;
}

void matrix_print(const Matrix* m, const char* matrix_label) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    fprintf(stderr, "==> ERRO: EMATRIXNULL");
    return;
  }

  unsigned int label_size = 0;

  if (matrix_label != NULL) {
    label_size = strlen(matrix_label) + 3;
    printf("%s = ", matrix_label);
  }

  for (unsigned int i = 0; i < m->nlins; i++) {
    if (i > 0 && label_size > 0) {
      for (unsigned int j = 0; j < label_size; j++) {
        printf(" ");
      }
    }

    printf("| ");

    for (unsigned int j = 0; j < m->ncols; j++) {
      printf("%.2f ", m->columns[j].x[i]);
    }

    printf("|\n");
  }
}

void matrix_println(const Matrix* m, const char* matrix_label) {
  matrix_print(m, matrix_label);
  printf("\n");
}
