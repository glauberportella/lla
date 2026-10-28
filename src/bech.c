/*
 * bech.c - Benchmark de performance da biblioteca LLA (Little Linear Algebra).
 *
 * Executa as operacoes de matrizes quadradas em dimensoes crescentes
 * (5, 10, 20 e 50 por padrao) e reporta tempo, throughput e MFLOP/s de cada
 * operacao.
 *
 * Uso:
 *   ./bech              dimensoes padrao: 5 10 20 50
 *   ./bech 3 7 100      dimensoes customizadas
 *
 * Observacoes:
 *   - Os tempos include o custo de alocacao e liberacao de memoria, ja que
 *     a biblioteca aloca o resultado de cada operacao.
 *   - As matrizes de entrada sao geradas de forma pseudoaleatoria e
 *     diagonalmente dominantes, o que garante que nao sejam singulares
 *     (necessario para determinante e inversa).
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "lu.h"
#include "matrix.h"
#include "vector.h"

#define BENCH_TARGET_SECONDS 0.20
#define BENCH_MAX_ITERATIONS 2000000u
#define BENCH_MIN_DIMENSION 1u
#define BENCH_MAX_DIMENSION 512u
#define BENCH_MAX_DIM_COUNT 16
#define BENCH_DEFAULT_DIM_COUNT 4

static const unsigned int bench_default_dims[BENCH_DEFAULT_DIM_COUNT] = { 5, 10, 20, 50 };

/**
 * Tempo monotono em segundos.
 *
 * Usa CLOCK_MONOTONIC quando disponivel e CLOCKS_PER_SEC como reserva.
 */
static double bench_now(void) {
#if defined(CLOCK_MONOTONIC)
  struct timespec ts;
  if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
  }
#endif
  return (double)clock() / (double)CLOCKS_PER_SEC;
}

/**
 * Gerador pseudoaleatorio xorshift64* com semente fixa, de forma que os
 * resultados do benchmark sejam reproduziveis entre execucoes.
 */
static unsigned long long bench_rng_state = 88172645463325252ULL;

static void bench_rng_seed(unsigned long long seed) {
  bench_rng_state = (seed != 0ULL) ? seed : 88172645463325252ULL;
}

static unsigned int bench_rng_next(void) {
  unsigned long long x = bench_rng_state;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  bench_rng_state = x;
  return (unsigned int)((x * 2685821657736338717ULL) >> 32);
}

/**
 * Valor double no intervalo [-1, 1).
 */
static double bench_rng_double(void) {
  return ((double)bench_rng_next() / 4294967296.0) * 2.0 - 1.0;
}

/**
 * Matriz quadrada n x n diagonalmente dominante, logo nao singular.
 */
static Matrix* bench_random_matrix(unsigned int n) {
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
      const double v = bench_rng_double();
      matrix_set(m, i, j, v);
      soma += fabs(v);
    }

    // A diagonal domina a linha, o que torna a matriz nao singular.
    matrix_set(m, i, i, soma + 1.0 + bench_rng_double());
  }

  return m;
}

/**
 * Soma dos elementos da matriz, lida direto do bloco denso para nao
 * distorcer a medicao de tempo com o custo da API matrix_get.
 */
static double bench_matrix_checksum(const Matrix* m) {
  if (m == NULL) {
    return NAN;
  }

  const size_t total = (size_t)m->nlins * (size_t)m->ncols;
  double acc = 0.0;
  for (size_t k = 0; k < total; k++) {
    acc += m->data[k];
  }

  return acc;
}

typedef double (*bench_op_fn)(const Matrix* a, const Matrix* b);
typedef double (*bench_flops_fn)(unsigned int n);

typedef struct {
  const char* name;
  bench_op_fn run;
  bench_flops_fn flops;
  int needs_b;
} BenchOp;

static double bench_op_mul(const Matrix* a, const Matrix* b) {
  Matrix* r = matrix_mul(a, b);
  if (r == NULL) {
    return NAN;
  }
  const double acc = bench_matrix_checksum(r);
  matrix_destroy(r);
  return acc;
}

static double bench_op_add(const Matrix* a, const Matrix* b) {
  Matrix* r = matrix_add(a, b);
  if (r == NULL) {
    return NAN;
  }
  const double acc = bench_matrix_checksum(r);
  matrix_destroy(r);
  return acc;
}

static double bench_op_sub(const Matrix* a, const Matrix* b) {
  Matrix* r = matrix_sub(a, b);
  if (r == NULL) {
    return NAN;
  }
  const double acc = bench_matrix_checksum(r);
  matrix_destroy(r);
  return acc;
}

static double bench_op_mul_scalar(const Matrix* a, const Matrix* b) {
  (void)b;
  Matrix* r = matrix_mul_scalar(2.5, a);
  if (r == NULL) {
    return NAN;
  }
  const double acc = bench_matrix_checksum(r);
  matrix_destroy(r);
  return acc;
}

static double bench_op_transpose(const Matrix* a, const Matrix* b) {
  (void)b;
  Matrix* r = matrix_transpose(a);
  if (r == NULL) {
    return NAN;
  }
  const double acc = bench_matrix_checksum(r);
  matrix_destroy(r);
  return acc;
}

static double bench_op_determinant(const Matrix* a, const Matrix* b) {
  (void)b;
  return matrix_determinant(a);
}

static double bench_op_inverse(const Matrix* a, const Matrix* b) {
  (void)b;
  Matrix* r = matrix_inverse(a);
  if (r == NULL) {
    return NAN;
  }
  const double acc = bench_matrix_checksum(r);
  matrix_destroy(r);
  return acc;
}

static double bench_op_lu(const Matrix* a, const Matrix* b) {
  (void)b;
  LU* lu = lu_create_decomposition(a);
  if (lu == NULL) {
    return NAN;
  }
  double acc = bench_matrix_checksum(lu->l);
  acc += bench_matrix_checksum(lu->u);
  acc += bench_matrix_checksum(lu->p);
  lu_destroy(lu);
  return acc;
}

static double bench_flops_mul(unsigned int n) {
  return 2.0 * (double)n * (double)n * (double)n;
}

static double bench_flops_quadratic(unsigned int n) {
  return (double)n * (double)n;
}

static double bench_flops_lu(unsigned int n) {
  return (2.0 / 3.0) * (double)n * (double)n * (double)n;
}

static double bench_flops_inverse(unsigned int n) {
  const double n3 = (double)n * (double)n * (double)n;
  return (2.0 / 3.0 + 2.0) * n3;
}

static const BenchOp bench_ops[] = {
  { "matrix_mul", bench_op_mul, bench_flops_mul, 1 },
  { "matrix_add", bench_op_add, bench_flops_quadratic, 1 },
  { "matrix_sub", bench_op_sub, bench_flops_quadratic, 1 },
  { "matrix_mul_scalar", bench_op_mul_scalar, bench_flops_quadratic, 0 },
  { "matrix_transpose", bench_op_transpose, bench_flops_quadratic, 0 },
  { "lu_decomposition", bench_op_lu, bench_flops_lu, 0 },
  { "matrix_determinant", bench_op_determinant, bench_flops_lu, 0 },
  { "matrix_inverse", bench_op_inverse, bench_flops_inverse, 0 },
};

#define BENCH_OP_COUNT (sizeof(bench_ops) / sizeof(bench_ops[0]))

typedef struct {
  unsigned int iterations;
  double total_seconds;
  double ms_per_op;
  double ops_per_second;
  double mflops;
  double checksum;
} BenchRow;

/**
 * Estima quantas repeticoes sao necessarias para atingir o tempo alvo,
 * medindo antes uma unica execucao da operacao.
 */
static unsigned int bench_calibrate(const BenchOp* op, const Matrix* a, const Matrix* b) {
  const double t0 = bench_now();
  const double warmup = op->run(a, b);
  const double t1 = bench_now();

  if (isnan(warmup)) {
    return 0;
  }

  double single = t1 - t0;
  if (single <= 0.0) {
    single = 1e-9;
  }

  double iterations = BENCH_TARGET_SECONDS / single;
  if (iterations < 1.0) {
    iterations = 1.0;
  }
  if (iterations > (double)BENCH_MAX_ITERATIONS) {
    iterations = (double)BENCH_MAX_ITERATIONS;
  }

  return (unsigned int)iterations;
}

/**
 * Executa a operacao repetidamente e preenche a linha de resultado.
 */
static int bench_run_op(const BenchOp* op, const Matrix* a, const Matrix* b, BenchRow* row) {
  const unsigned int iterations = bench_calibrate(op, a, b);
  if (iterations == 0) {
    return -1;
  }

  double checksum = 0.0;

  const double t0 = bench_now();
  for (unsigned int i = 0; i < iterations; i++) {
    const double valor = op->run(a, b);
    if (isnan(valor)) {
      return -1;
    }
    checksum += valor;
  }
  const double t1 = bench_now();

  row->iterations = iterations;
  row->total_seconds = t1 - t0;
  row->ms_per_op = (row->total_seconds * 1000.0) / (double)iterations;
  row->ops_per_second = (double)iterations / row->total_seconds;
  row->mflops = (op->flops((unsigned int)a->nlins) * (double)iterations) / row->total_seconds / 1e6;
  row->checksum = checksum;

  return 0;
}

static void bench_print_header(void) {
  printf("%-20s %12s %12s %14s %14s %10s\n",
         "operacao", "iteracoes", "tempo (s)", "ms/op", "ops/s", "MFLOP/s");
  printf("---------------------------------------------------------------------------------------------------\n");
}

static void bench_print_row(const BenchOp* op, const BenchRow* row) {
  printf("%-20s %12u %12.6f %14.6f %14.1f %10.3f\n",
         op->name, row->iterations, row->total_seconds, row->ms_per_op,
         row->ops_per_second, row->mflops);
}

int main(int argc, char* argv[]) {
  unsigned int dims[BENCH_MAX_DIM_COUNT];
  unsigned int dim_count = 0;

  if (argc > 1) {
    for (int i = 1; i < argc && dim_count < BENCH_MAX_DIM_COUNT; i++) {
      char* end = NULL;
      errno = 0;
      const long valor = strtol(argv[i], &end, 10);

      if (errno != 0 || end == argv[i] || *end != '\0' ||
          valor < (long)BENCH_MIN_DIMENSION || valor > (long)BENCH_MAX_DIMENSION) {
        fprintf(stderr, "ERRO: dimensao invalida '%s' (esperado um inteiro entre %u e %u)\n",
                argv[i], BENCH_MIN_DIMENSION, BENCH_MAX_DIMENSION);
        return 1;
      }

      dims[dim_count++] = (unsigned int)valor;
    }
  } else {
    for (unsigned int i = 0; i < BENCH_DEFAULT_DIM_COUNT; i++) {
      dims[dim_count++] = bench_default_dims[i];
    }
  }

  if (dim_count == 0) {
    fprintf(stderr, "ERRO: nenhuma dimensao informada\n");
    return 1;
  }

  printf("=== LLA - Benchmark de matrizes quadradas ===\n");
#if defined(__VERSION__)
  printf("compilador: %s\n", __VERSION__);
#endif
  printf("dimensoes: ");
  for (unsigned int i = 0; i < dim_count; i++) {
    printf("%u%s", dims[i], (i + 1 < dim_count) ? ", " : "\n");
  }
  printf("tempo alvo por operacao: %.2f s\n", BENCH_TARGET_SECONDS);
  printf("os tempos incluem o custo de alocacao de memoria da biblioteca\n\n");

  BenchRow results[BENCH_MAX_DIM_COUNT][BENCH_OP_COUNT];

  for (unsigned int d = 0; d < dim_count; d++) {
    const unsigned int n = dims[d];

    bench_rng_seed(0x9E3779B97F4A7C15ULL + n);

    Matrix* a = bench_random_matrix(n);
    if (a == NULL) {
      fprintf(stderr, "ERRO: falha ao criar a matriz %u x %u (errno = %d)\n", n, n, errno);
      return 1;
    }

    Matrix* b = bench_random_matrix(n);
    if (b == NULL) {
      fprintf(stderr, "ERRO: falha ao criar a matriz %u x %u (errno = %d)\n", n, n, errno);
      matrix_destroy(a);
      return 1;
    }

    printf(">>> n = %u\n", n);
    bench_print_header();

    for (size_t o = 0; o < BENCH_OP_COUNT; o++) {
      if (bench_run_op(&bench_ops[o], a, b, &results[d][o]) != 0) {
        fprintf(stderr, "ERRO: falha em %s para n = %u (errno = %d)\n", bench_ops[o].name, n, errno);
        matrix_destroy(b);
        matrix_destroy(a);
        return 1;
      }
      bench_print_row(&bench_ops[o], &results[d][o]);
    }

    printf("\n");

    matrix_destroy(b);
    matrix_destroy(a);
  }

  printf("=== Escala: ms por operacao ===\n");
  printf("%-20s", "operacao");
  for (unsigned int d = 0; d < dim_count; d++) {
    printf(" %12u", dims[d]);
  }
  printf("\n");
  printf("---------------------------------------------------------------------------------------------------\n");

  for (size_t o = 0; o < BENCH_OP_COUNT; o++) {
    printf("%-20s", bench_ops[o].name);
    for (unsigned int d = 0; d < dim_count; d++) {
      printf(" %12.6f", results[d][o].ms_per_op);
    }
    printf("\n");
  }

  printf("\n=== Checksums (verificacao de consistencia) ===\n");
  for (size_t o = 0; o < BENCH_OP_COUNT; o++) {
    printf("%-20s", bench_ops[o].name);
    for (unsigned int d = 0; d < dim_count; d++) {
      printf(" %18.6e", results[d][o].checksum);
    }
    printf("\n");
  }

  return 0;
}
