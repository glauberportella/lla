#include "vector.h"
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <errno.h>
#include <stddef.h>

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

void vector_print(const Vector* v) {
  if (v == NULL) {
    errno = EVECTORNULL;
    fprintf(stderr, "==>ERRO: EVECTORNULL");
    return;
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

void vector_println(const Vector* v) {
  vector_print(v);
  printf("\n");
}

