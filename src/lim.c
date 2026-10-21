/*
 * lim.c - Descobre os limites praticos de dimensao da biblioteca LLA.
 *
 * Para cada operacao de matriz, procura a maior dimensao quadrada que ainda
 * pode ser executada de forma confiavel, considerando dois limites:
 *
 *   - tempo: a operacao precisa caber no orcamento de tempo definido
 *     (padrao de 10 s por execucao);
 *   - memoria: as alocacoes da biblioteca precisam ter sucesso.
 *
 * A busca dobra a dimensao a cada passo (8, 16, 32, ...) e usa o expoente
 * de complexidade medido entre duas execucoes para prever o custo da
 * proxima dimensao, evitando travamentos em operacoes O(n^3).
 *
 * Ao final, o resultado de cada operacao no maior tamanho suportado e
 * conferido contra uma implementacao de referencia independente, para
 * garantir que "nao quebrou" tambem significa "calculou certo".
 *
 * Uso:
 *   ./lim            # orcamento padrao de 10 s
 *   ./lim 2.5        # orcamento de 2,5 s por execucao
 */

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/resource.h>

#if defined(__APPLE__)
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

#include "lu.h"
#include "matrix.h"
#include "vector.h"

#define LIM_START_DIM      8u
#define LIM_DEFAULT_BUDGET 10.0
#define LIM_MEMORY_FRACTION 0.60

typedef enum {
  LIM_OP_MUL,
  LIM_OP_ADD,
  LIM_OP_SUB,
  LIM_OP_SCALAR,
  LIM_OP_TRANSPOSE,
  LIM_OP_LU,
  LIM_OP_DET,
  LIM_OP_INV
} LimOpId;

typedef int (*lim_op_fn)(const Matrix* a, const Matrix* b);

typedef struct {
  const char* name;
  LimOpId id;
  lim_op_fn run;
  int needs_b;
  double exponent; /* ordem de grandeza esperada: 2 ou 3 */
} LimOp;

/* ------------------------------------------------------------------ */
/* Utilidades de tempo e memoria                                       */
/* ------------------------------------------------------------------ */

static double lim_now(void) {
#if defined(CLOCK_MONOTONIC)
  struct timespec ts;
  if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
  }
#endif
  return (double)clock() / (double)CLOCKS_PER_SEC;
}

static size_t lim_max_rss(void) {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) != 0) {
    return 0;
  }
#if defined(__APPLE__)
  return (size_t)usage.ru_maxrss; /* bytes */
#else
  return (size_t)usage.ru_maxrss * 1024u; /* kibibytes */
#endif
}

static double lim_physical_memory(void) {
#if defined(__APPLE__)
  uint64_t total = 0;
  size_t size = sizeof(total);
  int mib[2] = { CTL_HW, HW_MEMSIZE };
  if (sysctl(mib, 2, &total, &size, NULL, 0) == 0) {
    return (double)total;
  }
  return 0.0;
#else
  const long pages = sysconf(_SC_PHYS_PAGES);
  const long page_size = sysconf(_SC_PAGE_SIZE);
  if (pages > 0 && page_size > 0) {
    return (double)pages * (double)page_size;
  }
  return 0.0;
#endif
}

static void lim_print_bytes(double bytes) {
  const double gb = bytes / (1024.0 * 1024.0 * 1024.0);
  const double mb = bytes / (1024.0 * 1024.0);
  if (gb >= 1.0) {
    printf("%.2f GB", gb);
  } else {
    printf("%.2f MB", mb);
  }
}

/* ------------------------------------------------------------------ */
/* Geracao de matrizes                                                 */
/* ------------------------------------------------------------------ */

static unsigned long long lim_rng_state = 0x2545F4914F6CDD1DULL;

static void lim_rng_seed(unsigned long long seed) {
  lim_rng_state = (seed != 0ULL) ? seed : 0x2545F4914F6CDD1DULL;
}

static unsigned int lim_rng_next(void) {
  unsigned long long x = lim_rng_state;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  lim_rng_state = x;
  return (unsigned int)((x * 2685821657736338717ULL) >> 32);
}

static double lim_rng_double(void) {
  return ((double)lim_rng_next() / 4294967296.0) * 2.0 - 1.0;
}

/**
 * Matriz n x n diagonalmente dominante, portanto nao singular.
 */
static Matrix* lim_random_matrix(unsigned int n) {
  Matrix* m = matrix_create_empty(n, n);
  if (m == NULL) {
    return NULL;
  }

  for (unsigned int i = 0; i < n; i++) {
    double soma = 0.0;
    for (unsigned int j = 0; j < n; j++) {
      if (i == j) {
        continue;
      }
      const double v = lim_rng_double();
      matrix_set(m, i, j, v);
      soma += fabs(v);
    }
    matrix_set(m, i, i, soma + 1.0 + lim_rng_double());
  }

  return m;
}

/* ------------------------------------------------------------------ */
/* Implementacao de referencia (matriz densa row-major, sem validacao) */
/* ------------------------------------------------------------------ */

typedef struct {
  unsigned int n;
  double* data;
} RefMatrix;

static RefMatrix* ref_create(unsigned int n) {
  RefMatrix* m = (RefMatrix*)malloc(sizeof(RefMatrix));
  if (m == NULL) {
    return NULL;
  }
  m->n = n;
  m->data = (double*)calloc((size_t)n * n, sizeof(double));
  if (m->data == NULL) {
    free(m);
    return NULL;
  }
  return m;
}

static void ref_destroy(RefMatrix* m) {
  if (m != NULL) {
    free(m->data);
    free(m);
  }
}

static double ref_get(const RefMatrix* m, unsigned int i, unsigned int j) {
  return m->data[(size_t)i * m->n + j];
}

static void ref_set(RefMatrix* m, unsigned int i, unsigned int j, double v) {
  m->data[(size_t)i * m->n + j] = v;
}

/**
 * Converte uma matriz da biblioteca (vetores coluna) para a referencia.
 */
static RefMatrix* ref_from_matrix(const Matrix* m) {
  RefMatrix* r = ref_create(m->nlins);
  if (r == NULL) {
    return NULL;
  }
  for (unsigned int i = 0; i < m->nlins; i++) {
    for (unsigned int j = 0; j < m->ncols; j++) {
      ref_set(r, i, j, m->columns[j].x[i]);
    }
  }
  return r;
}

static RefMatrix* ref_mul(const RefMatrix* a, const RefMatrix* b) {
  if (a->n != b->n) {
    return NULL;
  }
  RefMatrix* c = ref_create(a->n);
  if (c == NULL) {
    return NULL;
  }
  for (unsigned int i = 0; i < a->n; i++) {
    for (unsigned int j = 0; j < b->n; j++) {
      double acc = 0.0;
      for (unsigned int k = 0; k < a->n; k++) {
        acc += ref_get(a, i, k) * ref_get(b, k, j);
      }
      ref_set(c, i, j, acc);
    }
  }
  return c;
}

/**
 * Determinante por eliminacao de Gauss com pivoteamento parcial.
 */
static double ref_determinant(const RefMatrix* m) {
  const unsigned int n = m->n;
  double* a = (double*)malloc((size_t)n * n * sizeof(double));
  if (a == NULL) {
    return NAN;
  }
  memcpy(a, m->data, (size_t)n * n * sizeof(double));

  double det = 1.0;
  for (unsigned int k = 0; k < n; k++) {
    unsigned int piv = k;
    double maximo = fabs(a[(size_t)k * n + k]);
    for (unsigned int i = k + 1; i < n; i++) {
      const double v = fabs(a[(size_t)i * n + k]);
      if (v > maximo) {
        maximo = v;
        piv = i;
      }
    }

    if (maximo == 0.0) {
      free(a);
      return 0.0;
    }

    if (piv != k) {
      for (unsigned int j = 0; j < n; j++) {
        const double tmp = a[(size_t)k * n + j];
        a[(size_t)k * n + j] = a[(size_t)piv * n + j];
        a[(size_t)piv * n + j] = tmp;
      }
      det = -det;
    }

    det *= a[(size_t)k * n + k];

    for (unsigned int i = k + 1; i < n; i++) {
      const double fator = a[(size_t)i * n + k] / a[(size_t)k * n + k];
      if (fator == 0.0) {
        continue;
      }
      for (unsigned int j = k; j < n; j++) {
        a[(size_t)i * n + j] -= fator * a[(size_t)k * n + j];
      }
    }
  }

  free(a);
  return det;
}

static double ref_max_diff(const RefMatrix* a, const RefMatrix* b) {
  double diff = 0.0;
  for (size_t i = 0; i < (size_t)a->n * a->n; i++) {
    const double d = fabs(a->data[i] - b->data[i]);
    if (d > diff) {
      diff = d;
    }
  }
  return diff;
}

static double ref_scale(const RefMatrix* a) {
  double maior = 0.0;
  for (size_t i = 0; i < (size_t)a->n * a->n; i++) {
    const double v = fabs(a->data[i]);
    if (v > maior) {
      maior = v;
    }
  }
  return maior;
}

/* ------------------------------------------------------------------ */
/* Operacoes da biblioteca sob teste                                   */
/* ------------------------------------------------------------------ */

static int lim_run_mul(const Matrix* a, const Matrix* b) {
  Matrix* r = matrix_mul(a, b);
  if (r == NULL) {
    return -1;
  }
  matrix_destroy(r);
  return 0;
}

static int lim_run_add(const Matrix* a, const Matrix* b) {
  Matrix* r = matrix_add(a, b);
  if (r == NULL) {
    return -1;
  }
  matrix_destroy(r);
  return 0;
}

static int lim_run_sub(const Matrix* a, const Matrix* b) {
  Matrix* r = matrix_sub(a, b);
  if (r == NULL) {
    return -1;
  }
  matrix_destroy(r);
  return 0;
}

static int lim_run_scalar(const Matrix* a, const Matrix* b) {
  (void)b;
  Matrix* r = matrix_mul_scalar(2.5, a);
  if (r == NULL) {
    return -1;
  }
  matrix_destroy(r);
  return 0;
}

static int lim_run_transpose(const Matrix* a, const Matrix* b) {
  (void)b;
  Matrix* r = matrix_transpose(a);
  if (r == NULL) {
    return -1;
  }
  matrix_destroy(r);
  return 0;
}

static int lim_run_lu(const Matrix* a, const Matrix* b) {
  (void)b;
  LU* lu = lu_create_decomposition(a);
  if (lu == NULL) {
    return -1;
  }
  lu_destroy(lu);
  return 0;
}

static int lim_run_det(const Matrix* a, const Matrix* b) {
  (void)b;
  const double det = matrix_determinant(a);
  if (isnan(det)) {
    return -1;
  }
  return 0;
}

static int lim_run_inv(const Matrix* a, const Matrix* b) {
  (void)b;
  Matrix* r = matrix_inverse(a);
  if (r == NULL) {
    return -1;
  }
  matrix_destroy(r);
  return 0;
}

static const LimOp lim_ops[] = {
  { "matrix_mul", LIM_OP_MUL, lim_run_mul, 1, 3.0 },
  { "matrix_add", LIM_OP_ADD, lim_run_add, 1, 2.0 },
  { "matrix_sub", LIM_OP_SUB, lim_run_sub, 1, 2.0 },
  { "matrix_mul_scalar", LIM_OP_SCALAR, lim_run_scalar, 0, 2.0 },
  { "matrix_transpose", LIM_OP_TRANSPOSE, lim_run_transpose, 0, 2.0 },
  { "lu_decomposition", LIM_OP_LU, lim_run_lu, 0, 3.0 },
  { "matrix_determinant", LIM_OP_DET, lim_run_det, 0, 3.0 },
  { "matrix_inverse", LIM_OP_INV, lim_run_inv, 0, 3.0 },
};

#define LIM_OP_COUNT (sizeof(lim_ops) / sizeof(lim_ops[0]))

/* ------------------------------------------------------------------ */
/* Verificacao numerica contra a referencia                            */
/* ------------------------------------------------------------------ */

/**
 * Executa a operacao em n e compara o resultado com a referencia.
 *
 * Devolve 0 quando os resultados conferem e -1 caso contrario.
 * Em `overflow` indica se o resultado saiu da faixa finita dos doubles.
 */
static int lim_verify(const LimOp* op, unsigned int n, double* erro, int* overflow) {
  *erro = 0.0;
  *overflow = 0;

  lim_rng_seed(0x9E3779B97F4A7C15ULL + n);

  Matrix* a = lim_random_matrix(n);
  if (a == NULL) {
    return -1;
  }
  Matrix* b = lim_random_matrix(n);
  if (b == NULL) {
    matrix_destroy(a);
    return -1;
  }

  RefMatrix* ra = ref_from_matrix(a);
  RefMatrix* rb = ref_from_matrix(b);
  int ok = -1;

  if (ra == NULL || rb == NULL) {
    goto cleanup;
  }

  switch (op->id) {
    case LIM_OP_MUL: {
      Matrix* r = matrix_mul(a, b);
      RefMatrix* rr = ref_mul(ra, rb);
      if (r != NULL && rr != NULL) {
        RefMatrix* rl = ref_from_matrix(r);
        if (rl != NULL) {
          const double escala = fmax(1.0, ref_scale(rr));
          *erro = ref_max_diff(rl, rr) / escala;
          ok = 0;
          ref_destroy(rl);
        }
        ref_destroy(rr);
      }
      if (r != NULL) {
        matrix_destroy(r);
      }
      break;
    }
    case LIM_OP_ADD:
    case LIM_OP_SUB: {
      Matrix* r = (op->id == LIM_OP_ADD) ? matrix_add(a, b) : matrix_sub(a, b);
      RefMatrix* rr = ref_create(n);
      if (r != NULL && rr != NULL) {
        RefMatrix* rl = ref_from_matrix(r);
        if (rl != NULL) {
          for (unsigned int i = 0; i < n; i++) {
            for (unsigned int j = 0; j < n; j++) {
              const double esperado = (op->id == LIM_OP_ADD)
                                          ? ref_get(ra, i, j) + ref_get(rb, i, j)
                                          : ref_get(ra, i, j) - ref_get(rb, i, j);
              ref_set(rr, i, j, esperado);
            }
          }
          const double escala = fmax(1.0, ref_scale(rr));
          *erro = ref_max_diff(rl, rr) / escala;
          ok = 0;
          ref_destroy(rl);
        }
        ref_destroy(rr);
      }
      if (r != NULL) {
        matrix_destroy(r);
      }
      break;
    }
    case LIM_OP_SCALAR: {
      Matrix* r = matrix_mul_scalar(2.5, a);
      RefMatrix* rr = ref_create(n);
      if (r != NULL && rr != NULL) {
        RefMatrix* rl = ref_from_matrix(r);
        if (rl != NULL) {
          for (unsigned int i = 0; i < n; i++) {
            for (unsigned int j = 0; j < n; j++) {
              ref_set(rr, i, j, 2.5 * ref_get(ra, i, j));
            }
          }
          const double escala = fmax(1.0, ref_scale(rr));
          *erro = ref_max_diff(rl, rr) / escala;
          ok = 0;
          ref_destroy(rl);
        }
        ref_destroy(rr);
      }
      if (r != NULL) {
        matrix_destroy(r);
      }
      break;
    }
    case LIM_OP_TRANSPOSE: {
      Matrix* r = matrix_transpose(a);
      RefMatrix* rr = ref_create(n);
      if (r != NULL && rr != NULL) {
        RefMatrix* rl = ref_from_matrix(r);
        if (rl != NULL) {
          for (unsigned int i = 0; i < n; i++) {
            for (unsigned int j = 0; j < n; j++) {
              ref_set(rr, i, j, ref_get(ra, j, i));
            }
          }
          *erro = ref_max_diff(rl, rr);
          ok = 0;
          ref_destroy(rl);
        }
        ref_destroy(rr);
      }
      if (r != NULL) {
        matrix_destroy(r);
      }
      break;
    }
    case LIM_OP_LU: {
      LU* lu = lu_create_decomposition(a);
      RefMatrix* rl = (lu != NULL) ? ref_from_matrix(lu->l) : NULL;
      RefMatrix* ru = (lu != NULL) ? ref_from_matrix(lu->u) : NULL;
      RefMatrix* rp = (lu != NULL) ? ref_from_matrix(lu->p) : NULL;
      if (rl != NULL && ru != NULL && rp != NULL) {
        RefMatrix* lu_prod = ref_mul(rl, ru);
        RefMatrix* pa = ref_mul(rp, ra);
        if (lu_prod != NULL && pa != NULL) {
          const double escala = fmax(1.0, ref_scale(pa));
          *erro = ref_max_diff(lu_prod, pa) / escala;
          ok = 0;
        }
        if (lu_prod != NULL) {
          ref_destroy(lu_prod);
        }
        if (pa != NULL) {
          ref_destroy(pa);
        }
      }
      if (rl != NULL) {
        ref_destroy(rl);
      }
      if (ru != NULL) {
        ref_destroy(ru);
      }
      if (rp != NULL) {
        ref_destroy(rp);
      }
      if (lu != NULL) {
        lu_destroy(lu);
      }
      break;
    }
    case LIM_OP_DET: {
      const double det = matrix_determinant(a);
      const double ref = ref_determinant(ra);
      if (!isnan(det) && !isnan(ref)) {
        if (isinf(det) || isinf(ref)) {
          *overflow = 1;
          *erro = (isinf(det) && isinf(ref) && (det > 0.0) == (ref > 0.0)) ? 0.0 : INFINITY;
        } else {
          const double escala = fmax(1.0, fabs(ref));
          *erro = fabs(det - ref) / escala;
        }
        ok = 0;
      }
      break;
    }
    case LIM_OP_INV: {
      Matrix* inv = matrix_inverse(a);
      RefMatrix* rinv = (inv != NULL) ? ref_from_matrix(inv) : NULL;
      if (rinv != NULL) {
        RefMatrix* produto = ref_mul(ra, rinv);
        if (produto != NULL) {
          double erro_max = 0.0;
          for (unsigned int i = 0; i < n; i++) {
            for (unsigned int j = 0; j < n; j++) {
              const double esperado = (i == j) ? 1.0 : 0.0;
              const double d = fabs(ref_get(produto, i, j) - esperado);
              if (d > erro_max) {
                erro_max = d;
              }
            }
          }
          *erro = erro_max;
          ok = 0;
          ref_destroy(produto);
        }
        ref_destroy(rinv);
      }
      if (inv != NULL) {
        matrix_destroy(inv);
      }
      break;
    }
    default:
      break;
  }

cleanup:
  if (ra != NULL) {
    ref_destroy(ra);
  }
  if (rb != NULL) {
    ref_destroy(rb);
  }
  matrix_destroy(b);
  matrix_destroy(a);

  return ok;
}

/* ------------------------------------------------------------------ */
/* Busca do limite por operacao                                        */
/* ------------------------------------------------------------------ */

typedef struct {
  const char* name;
  unsigned int max_n;
  double time_at_max;
  double predicted_next;
  size_t rss_at_max;
  int verified;
  double verify_error;
  int overflow;
  char reason[80];
} LimSummary;

/**
 * Procura a maior dimensao que a operacao suporta dentro do orcamento.
 */
static void lim_search(const LimOp* op, double budget, double memoria_limite, LimSummary* out) {
  unsigned int n = LIM_START_DIM;
  unsigned int n_anterior = 0;
  double t_anterior = 0.0;
  double expoente = op->exponent;

  out->name = op->name;
  out->max_n = 0;
  out->time_at_max = 0.0;
  out->predicted_next = 0.0;
  out->rss_at_max = 0;
  out->verified = 0;
  out->verify_error = 0.0;
  out->overflow = 0;
  snprintf(out->reason, sizeof(out->reason), "%s", "desconhecido");

  for (;;) {
    const double bytes_necessarios = (double)n * (double)n * (double)sizeof(double);
    if (bytes_necessarios * 8.0 > memoria_limite) {
      snprintf(out->reason, sizeof(out->reason), "memoria: %u x %u excede o limite", n, n);
      break;
    }

    lim_rng_seed(0x9E3779B97F4A7C15ULL + n);

    Matrix* a = lim_random_matrix(n);
    if (a == NULL) {
      snprintf(out->reason, sizeof(out->reason), "memoria: falha ao criar a matriz (errno = %d)", errno);
      break;
    }

    Matrix* b = NULL;
    if (op->needs_b) {
      b = lim_random_matrix(n);
      if (b == NULL) {
        matrix_destroy(a);
        snprintf(out->reason, sizeof(out->reason), "memoria: falha ao criar a segunda matriz (errno = %d)", errno);
        break;
      }
    }

    const double t0 = lim_now();
    const int rc = op->run(a, b);
    const double t1 = lim_now();
    const double decorrido = t1 - t0;

    matrix_destroy(b);
    matrix_destroy(a);

    if (rc != 0) {
      snprintf(out->reason, sizeof(out->reason), "erro da biblioteca (errno = %d)", errno);
      break;
    }

    out->max_n = n;
    out->time_at_max = decorrido;
    out->rss_at_max = lim_max_rss();

    if (decorrido > budget) {
      snprintf(out->reason, sizeof(out->reason), "tempo: excedeu o orcamento de %.1f s", budget);
      break;
    }

    if (n_anterior > 0 && t_anterior > 0.0) {
      const double razao_n = (double)n / (double)n_anterior;
      const double razao_t = decorrido / t_anterior;
      if (razao_n > 1.0 && razao_t > 0.0) {
        expoente = log(razao_t) / log(razao_n);
        if (expoente < 1.0) {
          expoente = 1.0;
        }
        if (expoente > 4.0) {
          expoente = 4.0;
        }
      }
    }

    const unsigned int proximo = n * 2u;
    const double razao = (double)proximo / (double)n;
    const double previsto = decorrido * pow(razao, expoente);
    out->predicted_next = previsto;

    if (previsto > budget) {
      snprintf(out->reason, sizeof(out->reason), "tempo: %u previsto em %.2f s", proximo, previsto);
      break;
    }

    n_anterior = n;
    t_anterior = decorrido;
    n = proximo;
  }

  if (out->max_n > 0) {
    double erro = 0.0;
    int overflow = 0;
    if (lim_verify(op, out->max_n, &erro, &overflow) == 0) {
      out->verified = 1;
      out->verify_error = erro;
      out->overflow = overflow;
    }
  }
}

/**
 * Avalia o determinante de uma matriz aleatoria de dimensao n.
 *
 * Devolve 1 quando o resultado e finito e 0 quando estourou a faixa.
 */
static int lim_det_is_finite(unsigned int n) {
  lim_rng_seed(0x9E3779B97F4A7C15ULL + n);

  Matrix* a = lim_random_matrix(n);
  if (a == NULL) {
    return 0;
  }

  const double det = matrix_determinant(a);
  matrix_destroy(a);

  return isfinite(det) ? 1 : 0;
}

/**
 * Procura a maior dimensao para a qual o determinante ainda cabe na faixa
 * finita dos doubles (o determinante cresce exponencialmente com n).
 *
 * Primeiro dobra a dimensao ate encontrar um determinante infinito e depois
 * refina o limiar com busca binaria.
 */
static unsigned int lim_det_finite_limit(double budget, double memoria_limite) {
  unsigned int n = LIM_START_DIM;
  unsigned int ultimo_finito = 0;
  unsigned int primeiro_infinito = 0;

  for (;;) {
    const unsigned int proximo = n * 2u;
    const double bytes = (double)proximo * (double)proximo * (double)sizeof(double) * 4.0;
    if (bytes > memoria_limite) {
      break;
    }

    if (!lim_det_is_finite(n)) {
      primeiro_infinito = n;
      break;
    }

    ultimo_finito = n;
    n = proximo;
  }

  if (primeiro_infinito == 0) {
    return ultimo_finito;
  }

  /* Refina o limiar: ultimo_finito e finito, primeiro_infinito nao e */
  unsigned int lo = ultimo_finito;
  unsigned int hi = primeiro_infinito;

  while (hi - lo > 1u) {
    const unsigned int mid = lo + (hi - lo) / 2u;
    const double bytes_mid = (double)mid * (double)mid * (double)sizeof(double) * 4.0;
    if (bytes_mid > memoria_limite) {
      break;
    }

    const double t0 = lim_now();
    const int finito = lim_det_is_finite(mid);
    const double t1 = lim_now();

    if (finito) {
      lo = mid;
    } else {
      hi = mid;
    }

    if ((t1 - t0) > budget) {
      break;
    }
  }

  return lo;
}

/**
 * Procura a maior dimensao que cabe na memoria, sem executar a operacao.
 */
static unsigned int lim_memory_search(unsigned int matrizes, double memoria_limite) {
  unsigned int n = LIM_START_DIM;

  for (;;) {
    const unsigned int proximo = n * 2u;
    const double bytes = (double)proximo * (double)proximo * (double)sizeof(double) * (double)matrizes;
    if (bytes > memoria_limite) {
      return n;
    }

    /* Confirma que a alocacao realmente funciona */
    Matrix** pilha = (Matrix**)calloc(matrizes, sizeof(Matrix*));
    if (pilha == NULL) {
      return n;
    }

    int ok = 1;
    for (unsigned int i = 0; i < matrizes; i++) {
      pilha[i] = matrix_create_empty(proximo, proximo);
      if (pilha[i] == NULL) {
        ok = 0;
        break;
      }
    }

    for (unsigned int i = 0; i < matrizes; i++) {
      matrix_destroy(pilha[i]);
    }
    free(pilha);

    if (!ok) {
      return n;
    }

    n = proximo;
  }
}

int main(int argc, char* argv[]) {
  double budget = LIM_DEFAULT_BUDGET;

  if (argc > 1) {
    char* end = NULL;
    errno = 0;
    const double valor = strtod(argv[1], &end);
    if (errno != 0 || end == argv[1] || *end != '\0' || valor <= 0.0) {
      fprintf(stderr, "ERRO: orcamento invalido '%s' (use um numero positivo de segundos)\n", argv[1]);
      return 1;
    }
    budget = valor;
  }

  const double fisica = lim_physical_memory();
  const double memoria_limite = (fisica > 0.0) ? fisica * LIM_MEMORY_FRACTION : 4.0 * 1024.0 * 1024.0 * 1024.0;

  printf("=== LLA - Limites de dimensao ===\n");
  if (fisica > 0.0) {
    printf("memoria fisica: ");
    lim_print_bytes(fisica);
    printf(" (usando %.0f%% como limite)\n", LIM_MEMORY_FRACTION * 100.0);
  }
  printf("orcamento de tempo por execucao: %.1f s\n", budget);
  printf("dimensao inicial: %u (dobrando a cada passo)\n\n", LIM_START_DIM);

  printf("=== Limite de memoria (somente alocacao) ===\n");
  const unsigned int n1 = lim_memory_search(1, memoria_limite);
  const unsigned int n3 = lim_memory_search(3, memoria_limite);
  const unsigned int n6 = lim_memory_search(6, memoria_limite);
  printf("1 matriz  n x n : %6u x %6u  (", n1, n1);
  lim_print_bytes((double)n1 * n1 * sizeof(double));
  printf(")\n");
  printf("3 matrizes      : %6u x %6u  (", n3, n3);
  lim_print_bytes((double)n3 * n3 * sizeof(double) * 3.0);
  printf(")\n");
  printf("6 matrizes      : %6u x %6u  (", n6, n6);
  lim_print_bytes((double)n6 * n6 * sizeof(double) * 6.0);
  printf(")\n");
  printf("3 matrizes e o pior caso do produto de matrizes;\n");
  printf("6 matrizes cobre a inversa (L, U, P, inversa, identidade e auxiliares).\n\n");

  LimSummary resumo[LIM_OP_COUNT];

  printf("=== Limite por operacao ===\n");
  for (size_t o = 0; o < LIM_OP_COUNT; o++) {
    printf("... %s\n", lim_ops[o].name);
    fflush(stdout);
    lim_search(&lim_ops[o], budget, memoria_limite, &resumo[o]);
  }

  printf("\n%-20s %10s %14s %12s %12s %s\n",
         "operacao", "max n", "tempo (s)", "previsto+", "pico RSS", "motivo");
  printf("---------------------------------------------------------------------------------------------------\n");
  for (size_t o = 0; o < LIM_OP_COUNT; o++) {
    printf("%-20s %10u %14.3f %12.2f ", resumo[o].name, resumo[o].max_n,
           resumo[o].time_at_max, resumo[o].predicted_next);
    lim_print_bytes((double)resumo[o].rss_at_max);
    printf(" %s\n", resumo[o].reason);
  }
  printf("(o pico RSS e a marca d'agua do processo, acumulada entre as operacoes)\n");

  printf("\n=== Verificacao numerica no maior tamanho suportado ===\n");
  printf("%-20s %10s %14s %s\n", "operacao", "n", "erro rel.", "resultado");
  printf("---------------------------------------------------------------------------------------------------\n");
  for (size_t o = 0; o < LIM_OP_COUNT; o++) {
    if (!resumo[o].verified) {
      printf("%-20s %10u %14s %s\n", resumo[o].name, resumo[o].max_n, "-", "FALHOU");
    } else if (resumo[o].overflow) {
      printf("%-20s %10u %14s %s\n", resumo[o].name, resumo[o].max_n, "-",
             "OVERFLOW (fora da faixa finita)");
    } else {
      printf("%-20s %10u %14.3e OK\n", resumo[o].name, resumo[o].max_n, resumo[o].verify_error);
    }
  }

  const unsigned int det_finito = lim_det_finite_limit(budget, memoria_limite);
  printf("\n=== Faixa finita do determinante ===\n");
  printf("o determinante cresce exponencialmente com n: para as matrizes usadas aqui\n");
  printf("ele deixa de caber em um double (|det| > 1.8e308) a partir de n = %u\n", det_finito + 1);
  printf("acima disso matrix_determinant devolve inf, e nao um numero valido\n");

  return 0;
}
