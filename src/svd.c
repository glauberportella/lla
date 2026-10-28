#include "svd.h"

#include <errno.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define SVD_TOLERANCE 1e-14

int lla_signal(double n) {
  return n >= 0 ? 1 : -1;
}

void svd_givens_rotation(double f, double g, double* c, double* s) {
  if (c == NULL || s == NULL) {
    errno = ESVDNULL;
    return;
  }

  const double norm = sqrt(f * f + g * g);
  if (norm == 0.0) {
    *c = 1.0;
    *s = 0.0;
    return;
  }

  *c = f / norm;
  *s = g / norm;
}

static double svd_dot(const double* a, const double* b, unsigned int n) {
  double sum = 0.0;
  for (unsigned int i = 0; i < n; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

/**
 * Monta o vetor de Householder v = x - alpha * e1 normalizado.
 *
 * O vetor e escrito em v (que pode ser o proprio x) e alpha (o novo valor
 * da diagonal) e devolvido. Retorna 0.0 quando o refletor e degenerado.
 */
static double svd_householder_vector(Vector* v, const double* x, unsigned int len) {
  const double norm_x = sqrt(svd_dot(x, x, len));
  if (norm_x == 0.0) {
    return 0.0;
  }

  const double alpha = -lla_signal(x[0]) * norm_x;

  for (unsigned int i = 0; i < len; i++) {
    v->x[i] = x[i];
  }
  v->x[0] -= alpha;

  const double norm_v = sqrt(svd_dot(v->x, v->x, len));
  if (norm_v == 0.0) {
    return 0.0;
  }

  for (unsigned int i = 0; i < len; i++) {
    v->x[i] /= norm_v;
  }

  return alpha;
}

/**
 * Aplica o refletor (I - 2 v v^T) as linhas k..m-1 de b e acumula o mesmo
 * refletor em u (u <- u * (I - 2 v v^T)).
 *
 * O produto v v^T nunca e materializado: cada coluna usa o produto interno
 * direto com v, o que torna a operacao O(m * n) e sem alocacoes.
 */
static void svd_apply_reflector_left(Matrix* b, Matrix* u, const Vector* v, unsigned int k) {
  const unsigned int m = b->nlins;
  const unsigned int n = b->ncols;

  for (unsigned int j = k; j < n; j++) {
    double* col = b->data + (size_t)j * b->nlins;
    const double dot = 2.0 * svd_dot(v->x, col + k, m - k);
    for (unsigned int i = k; i < m; i++) {
      col[i] -= v->x[i - k] * dot;
    }
  }

  const unsigned int ucols = u->ncols;
  for (unsigned int i = 0; i < m; i++) {
    double dot = 0.0;
    for (unsigned int l = k; l < ucols; l++) {
      dot += u->data[(size_t)l * u->nlins + i] * v->x[l - k];
    }
    dot *= 2.0;
    for (unsigned int j = k; j < ucols; j++) {
      u->data[(size_t)j * u->nlins + i] -= dot * v->x[j - k];
    }
  }
}

/**
 * Aplica o refletor (I - 2 v v^T) as colunas k+1..n-1 de b e acumula o
 * mesmo refletor em v_mat (v_mat <- v_mat * (I - 2 v v^T)).
 *
 * O produto v v^T nunca e materializado: cada linha usa o produto interno
 * direto com v, o que torna a operacao O(m * n) e sem alocacoes.
 */
static void svd_apply_reflector_right(Matrix* b, Matrix* v_mat, const Vector* v, unsigned int k) {
  const unsigned int m = b->nlins;
  const unsigned int n = b->ncols;

  for (unsigned int i = 0; i < m; i++) {
    double dot = 0.0;
    for (unsigned int l = k + 1; l < n; l++) {
      dot += b->data[(size_t)l * b->nlins + i] * v->x[l - k - 1];
    }
    dot *= 2.0;
    for (unsigned int j = k + 1; j < n; j++) {
      b->data[(size_t)j * b->nlins + i] -= dot * v->x[j - k - 1];
    }
  }

  for (unsigned int i = 0; i < v_mat->nlins; i++) {
    double dot = 0.0;
    for (unsigned int l = k + 1; l < n; l++) {
      dot += v->x[l - k - 1] * v_mat->data[(size_t)l * v_mat->nlins + i];
    }
    dot *= 2.0;
    for (unsigned int j = k + 1; j < n; j++) {
      v_mat->data[(size_t)j * v_mat->nlins + i] -= dot * v->x[j - k - 1];
    }
  }
}

/**
 * Rotacao a direita: v <- v * M, com M = [[c, s], [-s, c]] nas colunas p,p+1.
 */
static void svd_rotate_columns(Matrix* v, unsigned int p, double c, double s) {
  double* col_p = v->data + (size_t)p * v->nlins;
  double* col_p1 = v->data + (size_t)(p + 1) * v->nlins;

  for (unsigned int i = 0; i < v->nlins; i++) {
    const double vp = col_p[i];
    const double vp1 = col_p1[i];
    col_p[i] = c * vp - s * vp1;
    col_p1[i] = s * vp + c * vp1;
  }
}

/**
 * Rotacao a esquerda: u <- u * G, com G = [[c, s], [-s, c]] nas colunas p,p+1.
 */
static void svd_rotate_rows(Matrix* u, unsigned int p, double c, double s) {
  double* col_p = u->data + (size_t)p * u->nlins;
  double* col_p1 = u->data + (size_t)(p + 1) * u->nlins;

  for (unsigned int i = 0; i < u->nlins; i++) {
    const double up = col_p[i];
    const double up1 = col_p1[i];
    col_p[i] = c * up - s * up1;
    col_p1[i] = s * up + c * up1;
  }
}

static void svd_swap_columns(Matrix* m, unsigned int i, unsigned int j) {
  double* col_i = m->data + (size_t)i * m->nlins;
  double* col_j = m->data + (size_t)j * m->nlins;

  for (unsigned int l = 0; l < m->nlins; l++) {
    const double aux = col_i[l];
    col_i[l] = col_j[l];
    col_j[l] = aux;
  }
}

/**
 * Copia as ncols primeiras colunas de m.
 */
static Matrix* svd_trim_columns(const Matrix* m, unsigned int ncols) {
  Matrix* t = matrix_create_empty(m->nlins, ncols);
  if (t == NULL) {
    errno = ESVDALLOCATION;
    return NULL;
  }

  for (unsigned int j = 0; j < ncols; j++) {
    memcpy(t->data + (size_t)j * t->nlins, m->data + (size_t)j * m->nlins, m->nlins * sizeof(double));
  }

  return t;
}

HouseholderBidiagResult* svd_householder_bidiag(const Matrix* a) {
  if (a == NULL) {
    errno = ESVDNULL;
    return NULL;
  }

  const unsigned int m = a->nlins;
  const unsigned int n = a->ncols;

  if (m == 0 || n == 0) {
    errno = ESVDINVALIDDIM;
    return NULL;
  }

  
  Matrix* b = matrix_copy(a);
  if (b == NULL) {
    errno = ESVDALLOCATION;
    return NULL;
  }

  if (b == NULL) {
    errno = ESVDALLOCATION;
    return NULL;
  }

  Matrix* u = matrix_create_identity(m, m);
  if (u == NULL) {
    matrix_destroy(b);
    errno = ESVDALLOCATION;
    return NULL;
  }

  Matrix* v = matrix_create_identity(n, n);
  if (v == NULL) {
    matrix_destroy(u);
    matrix_destroy(b);
    errno = ESVDALLOCATION;
    return NULL;
  }

  Vector* d = vector_create_empty(n);
  if (d == NULL) {
    matrix_destroy(v);
    matrix_destroy(u);
    matrix_destroy(b);
    errno = ESVDALLOCATION;
    return NULL;
  }

  Vector* e = vector_create_empty(n > 1 ? n - 1 : 1);
  if (e == NULL) {
    vector_destroy(d);
    matrix_destroy(v);
    matrix_destroy(u);
    matrix_destroy(b);
    errno = ESVDALLOCATION;
    return NULL;
  }

  Vector* reflector = vector_create_empty(m > n ? m : n);
  if (reflector == NULL) {
    vector_destroy(e);
    vector_destroy(d);
    matrix_destroy(v);
    matrix_destroy(u);
    matrix_destroy(b);
    errno = ESVDALLOCATION;
    return NULL;
  }

  const unsigned int steps = (m < n) ? m : n;

  for (unsigned int k = 0; k < steps; k++) {
    if (k + 1 < m) {
      const double alpha = svd_householder_vector(reflector, b->data + (size_t)k * b->nlins + k, m - k);
      if (alpha != 0.0) {
        svd_apply_reflector_left(b, u, reflector, k);
        b->data[(size_t)k * b->nlins + k] = alpha;
      }
    }
    d->x[k] = b->data[(size_t)k * b->nlins + k];

    if (k + 1 < n) {
      for (unsigned int j = k + 1; j < n; j++) {
        reflector->x[j - k - 1] = b->data[(size_t)j * b->nlins + k];
      }
      const double alpha = svd_householder_vector(reflector, reflector->x, n - k - 1);
      if (alpha != 0.0) {
        svd_apply_reflector_right(b, v, reflector, k);
        b->data[(size_t)(k+1) * b->nlins + k] = alpha;
      }
    }

    if (k + 1 < n) {
      e->x[k] = b->data[(size_t)(k+1) * b->nlins + k];
    }
  }

  matrix_destroy(b);
  vector_destroy(reflector);

  Matrix* u_thin = svd_trim_columns(u, n);
  matrix_destroy(u);
  if (u_thin == NULL) {
    vector_destroy(e);
    vector_destroy(d);
    matrix_destroy(v);
    errno = ESVDALLOCATION;
    return NULL;
  }

  HouseholderBidiagResult* result = (HouseholderBidiagResult*)malloc(sizeof(HouseholderBidiagResult));
  if (result == NULL) {
    vector_destroy(e);
    vector_destroy(d);
    matrix_destroy(v);
    matrix_destroy(u_thin);
    errno = ESVDALLOCATION;
    return NULL;
  }

  result->d = d;
  result->e = e;
  result->U = u_thin;
  result->V = v;

  return result;
}

QRBidiagResult* svd_qr_bidiag_iter(HouseholderBidiagResult* hbdr) {
  if (hbdr == NULL || hbdr->d == NULL || hbdr->U == NULL || hbdr->V == NULL) {
    errno = ESVDNULL;
    return NULL;
  }

  const unsigned int n = hbdr->d->dim;
  Vector* alpha = hbdr->d;
  Vector* beta = hbdr->e;
  Matrix* u = hbdr->U;
  Matrix* v = hbdr->V;

  const unsigned int max_iter = 30 * n + 10;
  int converged = (n < 2);

  for (unsigned int sweep = 0; sweep < max_iter && !converged; sweep++) {
    unsigned int q = n - 1;
    while (q > 0 && fabs(beta->x[q - 1]) < SVD_TOLERANCE * (fabs(alpha->x[q - 1]) + fabs(alpha->x[q]))) {
      beta->x[q - 1] = 0.0;
      q--;
    }

    if (q == 0) {
      converged = 1;
      break;
    }

    const double f = alpha->x[q - 1] * alpha->x[q - 1] + (q >= 2 ? beta->x[q - 2] * beta->x[q - 2] : 0.0);
    const double g = alpha->x[q - 1] * beta->x[q - 1];
    const double h = alpha->x[q] * alpha->x[q] + beta->x[q - 1] * beta->x[q - 1];
    const double delta = (f - h) / 2.0;
    const double denom = delta + lla_signal(delta) * sqrt(delta * delta + g * g);
    const double mu = (denom != 0.0) ? h - (g * g) / denom : h;

    const double t00 = alpha->x[0] * alpha->x[0];
    const double t10 = alpha->x[0] * beta->x[0];
    const double r0 = sqrt((t00 - mu) * (t00 - mu) + t10 * t10);
    if (r0 == 0.0) {
      converged = 1;
      break;
    }

    const double c0 = (t00 - mu) / r0;
    const double s0 = -t10 / r0;
    const double a0 = alpha->x[0];
    const double b0 = beta->x[0];
    const double a1 = alpha->x[1];

    svd_rotate_columns(v, 0, c0, s0);
    alpha->x[0] = c0 * a0 - s0 * b0;
    beta->x[0] = s0 * a0 + c0 * b0;
    alpha->x[1] = c0 * a1;

    double gamma = -s0 * a1;

    for (unsigned int p = 0; p < q; p++) {
      const double ap = alpha->x[p];
      const double r = sqrt(ap * ap + gamma * gamma);
      double gbulge = 0.0;

      double c;
      double s;
      if (r > 0.0) {
        c = ap / r;
        s = -gamma / r;
      } else {
        c = 0.0;
        s = 1.0;
      }

      const double ap1 = alpha->x[p + 1];
      const double bp = beta->x[p];

      svd_rotate_rows(u, p, c, s);
      alpha->x[p] = r;
      beta->x[p] = c * bp - s * ap1;
      alpha->x[p + 1] = s * bp + c * ap1;

      if (p + 1 < n - 1) {
        const double bp1 = beta->x[p + 1];
        beta->x[p + 1] = c * bp1;
        gbulge = -s * bp1;
      }

      if (p + 2 <= q) {
        const double bp = beta->x[p];
        const double r2 = sqrt(bp * bp + gbulge * gbulge);

        if (r2 > 0.0) {
          const double c = bp / r2;
          const double s = -gbulge / r2;
          const double a1 = alpha->x[p + 1];
          const double b1 = beta->x[p + 1];
          const double a2 = alpha->x[p + 2];

          svd_rotate_columns(v, p + 1, c, s);
          beta->x[p] = r2;
          alpha->x[p + 1] = c * a1 - s * b1;
          beta->x[p + 1] = s * a1 + c * b1;
          alpha->x[p + 2] = c * a2;
          gamma = -s * a2;
        } else {
          gamma = 0.0;
        }
      }
    }

    for (unsigned int k = 0; k < q; k++) {
      if (fabs(beta->x[k]) < SVD_TOLERANCE * (fabs(alpha->x[k]) + fabs(alpha->x[k + 1]))) {
        beta->x[k] = 0.0;
      }
    }
  }

  if (!converged) {
    errno = ESVDNOTCONVERGED;
    return NULL;
  }

  QRBidiagResult* result = (QRBidiagResult*)malloc(sizeof(QRBidiagResult));
  if (result == NULL) {
    errno = ESVDALLOCATION;
    return NULL;
  }

  result->d = alpha;
  result->U = u;
  result->V = v;

  hbdr->d = NULL;
  hbdr->U = NULL;
  hbdr->V = NULL;
  vector_destroy(hbdr->e);
  hbdr->e = NULL;

  return result;
}

SVDResult* svd_pos_proc(QRBidiagResult* bdr) {
  if (bdr == NULL || bdr->d == NULL || bdr->U == NULL || bdr->V == NULL) {
    errno = ESVDNULL;
    return NULL;
  }

  const unsigned int n = bdr->d->dim;
  Vector* sigma = bdr->d;
  Matrix* u = bdr->U;
  Matrix* v = bdr->V;

  for (unsigned int i = 0; i < n; i++) {
    if (sigma->x[i] < 0.0) {
      sigma->x[i] = -sigma->x[i];
      for (unsigned int j = 0; j < u->nlins; j++) {
        u->data[(size_t)i * u->nlins + j] = -u->data[(size_t)i * u->nlins + j];
      }
    }
  }

  for (unsigned int i = 0; i + 1 < n; i++) {
    unsigned int maior = i;
    for (unsigned int j = i + 1; j < n; j++) {
      if (sigma->x[j] > sigma->x[maior]) {
        maior = j;
      }
    }

    if (maior != i) {
      const double aux = sigma->x[i];
      sigma->x[i] = sigma->x[maior];
      sigma->x[maior] = aux;
      svd_swap_columns(u, i, maior);
      svd_swap_columns(v, i, maior);
    }
  }

  Matrix* vt = matrix_transpose(v);
  if (vt == NULL) {
    errno = ESVDALLOCATION;
    return NULL;
  }
  matrix_destroy(v);

  SVDResult* result = (SVDResult*)malloc(sizeof(SVDResult));
  if (result == NULL) {
    matrix_destroy(vt);
    errno = ESVDALLOCATION;
    return NULL;
  }

  result->U = u;
  result->sigma = sigma;
  result->VT = vt;

  bdr->d = NULL;
  bdr->U = NULL;
  bdr->V = NULL;

  return result;
}

static SVDResult* svd_golub_reinsch_tall(const Matrix* a) {
  HouseholderBidiagResult* hbdr = svd_householder_bidiag(a);
  if (hbdr == NULL) {
    return NULL;
  }

  QRBidiagResult* bdr = svd_qr_bidiag_iter(hbdr);
  svd_householder_bidiag_destroy(hbdr);
  if (bdr == NULL) {
    return NULL;
  }

  SVDResult* result = svd_pos_proc(bdr);
  svd_qr_bidiag_destroy(bdr);
  return result;
}

SVDResult* svd_golub_reinsch(const Matrix* a) {
  if (a == NULL) {
    errno = ESVDNULL;
    return NULL;
  }

  const unsigned int m = a->nlins;
  const unsigned int n = a->ncols;

  if (m == 0 || n == 0) {
    errno = ESVDINVALIDDIM;
    return NULL;
  }

  if (m >= n) {
    return svd_golub_reinsch_tall(a);
  }

  Matrix* at = matrix_transpose(a);
  if (at == NULL) {
    errno = ESVDALLOCATION;
    return NULL;
  }

  SVDResult* at_result = svd_golub_reinsch_tall(at);
  matrix_destroy(at);
  if (at_result == NULL) {
    return NULL;
  }

  const unsigned int nsing = (m < n) ? m : n;

  Vector* sigma = vector_create_empty(nsing);
  if (sigma == NULL) {
    matrix_destroy(at_result->VT);
    matrix_destroy(at_result->U);
    vector_destroy(at_result->sigma);
    free(at_result);
    errno = ESVDALLOCATION;
    return NULL;
  }

  for (unsigned int i = 0; i < nsing; i++) {
    sigma->x[i] = at_result->sigma->x[i];
  }

  Matrix* u = matrix_transpose(at_result->VT);
  Matrix* vt = matrix_transpose(at_result->U);
  matrix_destroy(at_result->VT);
  matrix_destroy(at_result->U);
  vector_destroy(at_result->sigma);
  free(at_result);

  if (u == NULL || vt == NULL) {
    matrix_destroy(u);
    matrix_destroy(vt);
    vector_destroy(sigma);
    errno = ESVDALLOCATION;
    return NULL;
  }

  SVDResult* result = (SVDResult*)malloc(sizeof(SVDResult));
  if (result == NULL) {
    matrix_destroy(u);
    matrix_destroy(vt);
    vector_destroy(sigma);
    errno = ESVDALLOCATION;
    return NULL;
  }

  result->U = u;
  result->sigma = sigma;
  result->VT = vt;

  return result;
}

void svd_destroy(SVDResult* result) {
  if (result == NULL) {
    return;
  }

  matrix_destroy(result->U);
  vector_destroy(result->sigma);
  matrix_destroy(result->VT);
  free(result);
}

void svd_householder_bidiag_destroy(HouseholderBidiagResult* result) {
  if (result == NULL) {
    return;
  }

  vector_destroy(result->d);
  vector_destroy(result->e);
  matrix_destroy(result->U);
  matrix_destroy(result->V);
  free(result);
}

void svd_qr_bidiag_destroy(QRBidiagResult* result) {
  if (result == NULL) {
    return;
  }

  vector_destroy(result->d);
  matrix_destroy(result->U);
  matrix_destroy(result->V);
  free(result);
}
