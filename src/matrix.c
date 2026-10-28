#include "matrix.h"
#include "lu.h"
#include "svd.h"

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <float.h>

/*
 * Convenções internas do layout denso column-major:
 *
 *   a(i, j) = m->data[(size_t)j * m->nlins + i]
 *
 * Ou seja, cada coluna j é o trecho contíguo
 *
 *   m->data + j * nlins   ...   m->data + j * nlins + nlins - 1
 *
 * e a linha i tem passo (stride) nlins. O acesso às colunas é contíguo,
 * o que aproveita bem a cache e permite ao compilador vetorizar os
 * laços internos das operações.
 */

static size_t matrix_data_size(const Matrix* m) {
  return (size_t)m->nlins * (size_t)m->ncols;
}

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
  m->data = (double*)calloc(matrix_data_size(m), sizeof(double));

  if (m->data == NULL) {
    free(m);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  return m;
}

Matrix* matrix_create(unsigned int nlins, unsigned int ncols, Vector columns[]) {
  if (nlins == 0 || ncols == 0 || columns == NULL) {
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
  m->data = (double*)malloc(matrix_data_size(m) * sizeof(double));

  if (m->data == NULL) {
    free(m);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  for (unsigned int j = 0; j < ncols; j++) {
    if (columns[j].dim != nlins || columns[j].x == NULL) {
      free(m->data);
      free(m);
      errno = EMATRIXALLOCATION;
      return NULL;
    }

    memcpy(m->data + (size_t)j * nlins, columns[j].x, nlins * sizeof(double));
  }

  return m;
}

Matrix* matrix_create_identity(unsigned int m, unsigned int n) {
  Matrix* a = matrix_create_empty(m, n);
  if (a == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  const unsigned int diag = (m < n) ? m : n;
  for (unsigned int i = 0; i < diag; i++) {
    a->data[(size_t)i * n + i] = 1.0;
  }

  return a;
}

Matrix* matrix_copy(const Matrix* m) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  Matrix* c = matrix_create_empty(m->nlins, m->ncols);
  if (c == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  memcpy(c->data, m->data, matrix_data_size(m) * sizeof(double));

  return c;
}

void matrix_destroy(Matrix* m) {
  if (m != NULL) {
    free(m->data);
    free(m);
  }
}

int matrix_validate_index(const Matrix* m, const unsigned int i, const unsigned int j) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return -2;
  }

  if (i >= m->nlins) {
    errno = EMATRIXINVALIDLINEINDEX;
    return -1;
  }

  if (j >= m->ncols) {
    errno = EMATRIXINVALIDCOLINDEX;
    return 1;
  }

  return 0;
}

Vector* matrix_get_line_vectors(const Matrix* a) {
  if (a == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  Vector* line_vectors = (Vector*)calloc(a->nlins, sizeof(Vector));
  if (line_vectors == NULL) {
    errno = EVECTORALLOCATION;
    return NULL;
  }

  for (unsigned int j = 0; j < a->nlins; j++) {
    line_vectors[j].dim = a->ncols;
    line_vectors[j].x = (double*)calloc(a->ncols, sizeof(double));

    if (line_vectors[j].x == NULL) {
      for (unsigned int k = 0; k < j; k++) {
        free(line_vectors[k].x);
      }

      free(line_vectors);
      errno = EVECTORALLOCATION;
      return NULL;
    }

    for (unsigned int i = 0; i < a->ncols; i++) {
      line_vectors[j].x[i] = a->data[(size_t)i * a->nlins + j];
    }
  }

  return line_vectors;
}

Vector* matrix_get_vector_column_at(const Matrix* a, unsigned int i) {
  if (a == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  if (i >= a->ncols) {
    errno = EMATRIXINVALIDCOLINDEX;
    return NULL;
  }

  return vector_create(a->nlins, a->data + (size_t)i * a->nlins);
}

Matrix* matrix_set(Matrix* m, const unsigned int i, const unsigned int j, const double aij) {
  if (m == NULL || m->data == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  const int validated = matrix_validate_index(m, i, j);
  if (validated != 0) {
    return NULL;
  }

  m->data[(size_t)j * m->nlins + i] = aij;

  return m;
}

double matrix_get(const Matrix* m, const unsigned int i, const unsigned int j) {
  if (m == NULL || m->data == NULL) {
    errno = EMATRIXNULL;
    return NAN;
  }

  const int validated = matrix_validate_index(m, i, j);
  if (validated != 0) {
    return NAN;
  }

  return m->data[(size_t)j * m->nlins + i];
}

Matrix* matrix_remove_column_vector_at(const Matrix* m, unsigned int column_index) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  if (column_index > m->ncols - 1) {
    errno = EMATRIXINVALIDCOLINDEX;
    return NULL;
  }

  Matrix* n = matrix_create_empty(m->nlins, m->ncols - 1);
  if (n == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  const size_t col_bytes = m->nlins * sizeof(double);
  double* dest = n->data;
  for (unsigned int i = 0; i < m->ncols; i++) {
    if (i == column_index) {
      continue;
    }
    memcpy(dest, m->data + (size_t)i * m->nlins, col_bytes);
    dest += m->nlins;
  }

  return n;
}

/**
 * Substitui os valores da coluna de indice i pelos xlen valores de x.
 *
 * Os valores sao COPIADOS para dentro do bloco da matriz; o chamador
 * continua responsavel por liberar x.
 */
void matrix_add_column_vector_at(Matrix* m, unsigned int i, unsigned int xlen, double* x) {
  if (m == NULL || m->data == NULL) {
    errno = EMATRIXNULL;
    return;
  }

  if (i >= m->ncols) {
    errno = EMATRIXINVALIDCOLINDEX;
    return;
  }

  if (x == NULL) {
    errno = EVECTORNULL;
    return;
  }

  if (xlen != m->nlins) {
    errno = EINVALIDDIM;
    return;
  }

  memcpy(m->data + (size_t)i * m->nlins, x, xlen * sizeof(double));
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

  Matrix* c = matrix_create_empty(a->nlins, a->ncols);
  if (c == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  const size_t total = matrix_data_size(a);
  const double* restrict pa = a->data;
  const double* restrict pb = b->data;
  double* restrict pc = c->data;

  for (size_t k = 0; k < total; k++) {
    pc[k] = pa[k] + pb[k];
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

  Matrix* c = matrix_create_empty(a->nlins, a->ncols);
  if (c == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  const size_t total = matrix_data_size(a);
  const double* restrict pa = a->data;
  const double* restrict pb = b->data;
  double* restrict pc = c->data;

  for (size_t k = 0; k < total; k++) {
    pc[k] = pa[k] - pb[k];
  }

  return c;
}

/*
 * Produto C = A * B em laço j -> k -> i (saxpy por colunas de C):
 *
 *   para cada coluna j de C:
 *     C_col_j = 0
 *     para k: C_col_j += A_col_k * b(k, j)
 *
 * Com column-major, C_col_j, A_col_k e o elemento b(k, j) são acessados
 * de forma contígua, e o laço interno (sobre i) é vetorizado pelo
 * compilador. Se três streams simultâneos não bastarem, o OpenMP é
 * compilado no paralelismo sobre j quando o build usa OPENMP=1.
 */
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

  const unsigned int M = a->nlins; // linhas de A e de C
  const unsigned int K = a->ncols; // colunas de A == linhas de B
  const unsigned int N = b->ncols; // colunas de B e de C

  const double* restrict pa = a->data;
  const double* restrict pb = b->data;
  double* restrict pc = c->data;

#ifdef LLA_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int j = 0; j < N; j++) {
    double* restrict cj = pc + (size_t)j * M;
    const double* restrict bj = pb + (size_t)j * K; // coluna j de B, contígua
    for (unsigned int i = 0; i < M; i++) {
      cj[i] = 0.0;
    }
    for (unsigned int k = 0; k < K; k++) {
      const double bkj = bj[k];
      const double* restrict ak = pa + (size_t)k * M;
      for (unsigned int i = 0; i < M; i++) {
        cj[i] += ak[i] * bkj;
      }
    }
  }

  return c;
}

Matrix* matrix_mul_scalar(const double scalar, const Matrix* m) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  Matrix* c = matrix_create_empty(m->nlins, m->ncols);
  if (c == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  const size_t total = matrix_data_size(m);
  const double* restrict pm = m->data;
  double* restrict pc = c->data;

  for (size_t k = 0; k < total; k++) {
    pc[k] = scalar * pm[k];
  }

  return c;
}

Vector* matrix_mul_vector(const Matrix* m, const Vector* v) {
  if (m == NULL || m->data == NULL || v == NULL || v->x == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  if (m->ncols != v->dim) {
    errno = EMATRIXINVALIDMUL;
    return NULL;
  }

  Vector* r = vector_create_empty(m->nlins);
  if (r == NULL) {
    errno = EVECTORALLOCATION;
    return NULL;
  }

  // r = soma_j v[j] * col_j, com cada coluna lida de forma contígua.
  const unsigned int nr = m->nlins;
  double* restrict pr = r->x;

  for (unsigned int j = 0; j < m->ncols; j++) {
    const double coef = v->x[j];
    if (coef == 0.0) {
      continue;
    }
    const double* restrict col = m->data + (size_t)j * nr;
    for (unsigned int i = 0; i < nr; i++) {
      pr[i] += coef * col[i];
    }
  }

  return r;
}

Matrix* matrix_transpose(const Matrix* m) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  Matrix* t = matrix_create_empty(m->ncols, m->nlins);
  if (t == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  // t(j, i) = m(i, j): lê colunas de m (contíguo) e escreve as linhas de t.
  const unsigned int nr = m->nlins;
  const unsigned int nc = m->ncols;
  const double* restrict pm = m->data;
  double* restrict pt = t->data;

  for (unsigned int i = 0; i < nc; i++) {
    const double* restrict col = pm + (size_t)i * nr;
    for (unsigned int j = 0; j < nr; j++) {
      pt[(size_t)j * nc + i] = col[j];
    }
  }

  return t;
}

/**
 * Substituição direta in-place: resolve L * y = pb, com L triangular
 * inferior de diagonal unitária, atualizando o vetor y (saxpy por colunas,
 * com as colunas de L lidas de forma contígua).
 */
static void matrix_lu_forward_update(double* restrict y, const Matrix* restrict l, unsigned int n) {
  for (unsigned int j = 0; j < n; j++) {
    const double yj = y[j];
    if (yj == 0.0) {
      continue;
    }
    const double* restrict col = l->data + (size_t)j * n;
    for (unsigned int i = j + 1; i < n; i++) {
      y[i] -= col[i] * yj;
    }
  }
}

/**
 * Substituição retroativa: resolve U * x = y, com U triangular superior,
 * usando a forma por colunas (x[j] = y[j] / u_jj seguida da correção de
 * y[i], i < j, com a coluna j de U lida de forma contígua).
 */
static int matrix_lu_back_solve(double* restrict x, double* restrict y, const Matrix* restrict u, unsigned int n) {
  for (int jj = (int)n - 1; jj >= 0; jj--) {
    const unsigned int j = (unsigned int)jj;
    const double pivot = u->data[(size_t)j * n + j];

    if (fabs(pivot) < 1e-9) {
      errno = EMATRIXSINGULAR;
      return -1;
    }

    const double xj = y[j] / pivot;
    x[j] = xj;

    if (xj == 0.0) {
      continue;
    }

    const double* restrict col = u->data + (size_t)j * n;
    for (unsigned int i = 0; i < j; i++) {
      y[i] -= col[i] * xj;
    }
  }

  return 0;
}

Matrix* matrix_inverse(const Matrix* m) {
  /*
  Entrada: Matrizes L, U, P de tamanho n x n (onde PA = LU)
  Saída: Matriz A_inv (inversa de A)

  1. Decompor via LU: P * A = L * U (o vetor de pivôs registra a permutação).
  2. Para cada coluna i da inversa:
     a. montar b = P * e_i a partir do vetor de pivôs (b[k] = 1 se e só se
        a linha k de A foi movida para a posição k de P; nenhum produto
        matricial é preciso);
     b. [Substituição Direta] Resolver L * y = b por colunas de L;
     c. [Substituição Retroativa] Resolver U * x = y por colunas de U;
     d. copiar x como a i-ésima coluna de A_inv (A_inv[:, i] = x).
  3. Retornar A_inv.

  Trata-se da mesma matemática da versão anterior, mas sem construir a
  matriz identidade, sem criar matrizes n x 1 por coluna e sem multiplicações
  por P: o custo dominante são as n substituições O(n^2), que é o custo
  teórico mínimo da inversa.
  */

  if (m == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  if (m->nlins != m->ncols) {
    errno = EMATRIXINVALIDDIM;
    return NULL;
  }

  // Decomposição LU
  LU* lu = lu_create_decomposition(m);
  if (lu == NULL) {
    // lu_create_decomposition já define o errno apropriado
    return NULL;
  }

  unsigned int n = m->nlins;

  // A inversa a ser retornada
  Matrix* inv = matrix_create_empty(n, n);
  if (inv == NULL) {
    lu_destroy(lu);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  double* y = (double*)malloc(n * sizeof(double));
  if (y == NULL) {
    matrix_destroy(inv);
    lu_destroy(lu);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  const Matrix* l = lu->l;
  const Matrix* u = lu->u;
  const int* restrict piv = lu->piv;
  int failed = 0;

  for (unsigned int i = 0; i < n && !failed; i++) {
    // b = P * e_i: 1 na linha cujo pivô aponta para a coluna i
    for (unsigned int k = 0; k < n; k++) {
      y[k] = (piv[k] == (int)i) ? 1.0 : 0.0;
    }

    // L * y = b, com L de diagonal unitária (só as atualizações)
    matrix_lu_forward_update(y, l, n);

    // U * x = y
    double* x = inv->data + (size_t)i * n;
    if (matrix_lu_back_solve(x, y, u, n) != 0) {
      failed = 1;
    }
  }

  free(y);
  lu_destroy(lu);

  if (failed) {
    matrix_destroy(inv);
    return NULL;
  }

  return inv;
}

/**
 * Calcula o sinal da permutacao registrada no vetor de pivots do LU.
 *
 * Retorna 1 para permutacao par e -1 para permutacao impar.
 */
static int matrix_permutation_sign(const int* restrict piv, const unsigned int n) {
  // O sinal da permutacao e (-1) elevado ao numero de inversoes
  unsigned long inversions = 0;
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = i + 1; j < n; j++) {
      if (piv[i] > piv[j]) {
        inversions++;
      }
    }
  }

  return (inversions % 2 == 0) ? 1 : -1;
}

double matrix_determinant(const Matrix* m) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return NAN;
  }

  if (m->nlins != m->ncols) {
    errno = EMATRIXINVALIDDIM;
    return NAN;
  }

  LU* lu = lu_create_decomposition(m);
  if (lu == NULL) {
    // lu_create_decomposition ja define o errno (EMATRIXSINGULAR, EMATRIXINVALIDDIM, ...)
    if (errno == EMATRIXSINGULAR) {
      // Matriz singular tem determinante zero
      return 0.0;
    }

    errno = EMATRIXINVALIDDETERMINANT;
    return NAN;
  }

  const unsigned int n = m->nlins;

  /*
   * Como P * A = L * U, temos A = P^-1 * L * U e portanto
   * det(A) = det(P) * det(L) * det(U) = det(P) * prod(U[i][i]),
   * ja que L tem a diagonal unitaria (det(L) = 1).
   */
  double product = 1.0;
  const double* restrict du = lu->u->data;
  for (unsigned int i = 0; i < n; i++) {
    product *= du[(size_t)i * n + i];
  }

  const int sign = matrix_permutation_sign(lu->piv, n);

  lu_destroy(lu);

  return sign * product;
}

double matrix_trace(const Matrix* m) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return NAN;
  }

  if (m->nlins != m->ncols) {
    errno = EMATRIXINVALIDDIM;
    return NAN;
  }

  double trace = 0.0;
  for (unsigned int i = 0; i < m->nlins; i++) {
    trace += m->data[(size_t)i * m->nlins + i];
  }

  return trace;
}

unsigned int matrix_rank(const Matrix* m) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return (unsigned int)-1;
  }

  if (m->nlins == 0 || m->ncols == 0) {
    errno = EMATRIXINVALIDDIM;
    return (unsigned int)-1;
  }

  SVDResult* r = svd_golub_reinsch(m);
  if (r == NULL) {
    return (unsigned int)-1;
  }

  const double tol = (double)(m->nlins > m->ncols ? m->nlins : m->ncols)
                     * r->sigma->x[0] * DBL_EPSILON;

  unsigned int rank = 0;
  for (unsigned int j = 0; j < r->sigma->dim; j++) {
    if (r->sigma->x[j] > tol) {
      rank++;
    }
  }

  svd_destroy(r);
  return rank;
}

void matrix_print(const Matrix* m, const char* matrix_label) {
  if (m == NULL || m->data == NULL) {
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
      printf("%.2f ", m->data[(size_t)j * m->nlins + i]);
    }

    printf("|\n");
  }
}

void matrix_println(const Matrix* m, const char* matrix_label) {
  matrix_print(m, matrix_label);
  printf("\n");
}
