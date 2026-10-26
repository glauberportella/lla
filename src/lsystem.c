#include "lsystem.h"
#include "lu.h"
#include "svd.h"

#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <errno.h>

/**
 * Resolve L * y = pb, com L triangular inferior de diagonal unitaria.
 */
static Vector* lsystem_forward_substitution(const Matrix* l, const Vector* pb) {
  Vector* y = vector_create_empty(pb->dim);
  if (y == NULL) {
    return NULL;
  }

  for (unsigned int i = 0; i < pb->dim; i++) {
    double s = 0.0;
    for (unsigned int j = 0; j < i; j++) {
      s += matrix_get(l, i, j) * y->x[j];
    }
    vector_set(y, i, pb->x[i] - s);
  }

  return y;
}

/**
 * Resolve U * x = y, com U triangular superior.
 */
static Vector* lsystem_back_substitution(const Matrix* u, const Vector* y) {
  Vector* x = vector_create_empty(y->dim);
  if (x == NULL) {
    return NULL;
  }

  for (unsigned int i = y->dim; i > 0; i--) {
    const unsigned int ii = i - 1;

    double s = 0.0;
    for (unsigned int j = ii + 1; j < y->dim; j++) {
      s += matrix_get(u, ii, j) * x->x[j];
    }

    const double pivot = matrix_get(u, ii, ii);
    if (pivot == 0.0) {
      vector_destroy(x);
      return NULL;
    }
    vector_set(x, ii, (y->x[ii] - s) / pivot);
  }

  return x;
}

/**
 * Resolve o sistema quadrado via decomposicao LU: P * A = L * U, e
 * A * x = b equivale a L * U * x = P * b.
 *
 * Retorna NULL com errno = ELSYSTEMINVALIDDIM se o triangulamento encontrar
 * pivots nulos (matriz singular). Caso o pivot seja apenas pequeno, ainda
 * triangula: a decisao de cair para a SVD e do chamador.
 */
static Vector* lsystem_solve_lu(const Matrix* a, const Vector* b) {
  LU* lu = lu_create_decomposition(a);
  if (lu == NULL) {
    return NULL;
  }

  Vector* pb = matrix_mul_vector(lu->p, b);
  if (pb == NULL) {
    lu_destroy(lu);
    return NULL;
  }

  Vector* y = lsystem_forward_substitution(lu->l, pb);
  vector_destroy(pb);
  if (y == NULL) {
    lu_destroy(lu);
    return NULL;
  }

  Vector* x = lsystem_back_substitution(lu->u, y);
  vector_destroy(y);
  lu_destroy(lu);

  if (x == NULL && errno == 0) {
    errno = ELSYSTEMINVALIDDIM;
  }
  return x;
}

/**
 * Solucao de norma minima dos minimos quadrados: x = V * Sigma+ * U^T * b.
 *
 * Valores singulares abaixo de eps * max(m, n) * sigma_max sao descartados.
 */
static Vector* lsystem_solve_svd(const Matrix* a, const Vector* b) {
  SVDResult* r = svd_golub_reinsch(a);
  if (r == NULL) {
    return NULL;
  }

  const unsigned int k = r->sigma->dim; // min(m, n)
  const unsigned int nc = a->ncols;     // dimensao da solucao

  double utb[k];
  for (unsigned int j = 0; j < k; j++) {
    double s = 0.0;
    for (unsigned int i = 0; i < b->dim; i++) {
      s += matrix_get(r->U, i, j) * b->x[i];
    }
    utb[j] = s;
  }

  const double tol = (double)(a->nlins > nc ? a->nlins : nc)
                     * r->sigma->x[0] * DBL_EPSILON;

  Vector* x = vector_create_empty(nc);
  if (x == NULL) {
    svd_destroy(r);
    return NULL;
  }

  for (unsigned int j = 0; j < k; j++) {
    if (r->sigma->x[j] <= tol) {
      continue;
    }
    const double w = utb[j] / r->sigma->x[j];
    for (unsigned int i = 0; i < nc; i++) {
      x->x[i] += matrix_get(r->VT, j, i) * w;
    }
  }

  svd_destroy(r);
  return x;
}

Vector* lsystem_solve(const Matrix* a, const Vector* b) {
  if (a == NULL || b == NULL) {
    errno = ELSYSTEMNULL;
    return NULL;
  }

  if (a->nlins == 0 || a->ncols == 0 || b->dim != a->nlins) {
    errno = ELSYSTEMINVALIDDIM;
    return NULL;
  }

  if (a->nlins == a->ncols) {
    Vector* x = lsystem_solve_lu(a, b);
    if (x != NULL) {
      return x;
    }
    /* triangulamento com pivot nulo: cai para a SVD */
  }

  return lsystem_solve_svd(a, b);
}
