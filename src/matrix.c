#include "matrix.h"
#include "lu.h"

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
        free(m->columns[j].x);
      }
      free(m->columns);
      free(m);
      errno = EMATRIXALLOCATION;
      return NULL;
    }

    m->columns[i].dim = columns[i].dim;
    m->columns[i].x = (double*)malloc(m->columns[i].dim * sizeof(double));
    if (m->columns[i].x == NULL) {
      for (unsigned int j = 0; j < i; j++) {
        free(m->columns[j].x);
      }
      free(m->columns);
      free(m);
      errno = EMATRIXALLOCATION;
      return NULL;
    }

    memcpy(m->columns[i].x, columns[i].x, columns[i].dim * sizeof(double));
  }

  return m;
}

Matrix* matrix_create_identity(unsigned int m, unsigned int n) {
  Matrix* a = matrix_create_empty(m, n);
  if (a == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }
  for (unsigned int i = 0; i < m; i++) {
    for (unsigned int j = 0; j < n; j++) {
      if (i == j) {
        matrix_set(a, i, j, 1.0);
      }
    }
  }
  return a;
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

  if (a->columns == NULL) {
    errno = EVECTORNULL;
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
      line_vectors[j].x[i] = a->columns[i].x[j];
    }
  }

  return line_vectors;
}

Vector* matrix_get_vector_column_at(const Matrix* a, unsigned int i) {
  if (a == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }
  return &a->columns[i];
}

/**
 * Libera um conjunto de vetores linha devolvido por matrix_get_line_vectors.
 */
static void matrix_free_line_vectors(Vector* line_vectors, const unsigned int n) {
  if (line_vectors == NULL) {
    return;
  }

  for (unsigned int i = 0; i < n; i++) {
    if (line_vectors[i].x != NULL) {
      free(line_vectors[i].x);
    }
  }

  free(line_vectors);
}

Matrix* matrix_set(Matrix* m, const unsigned int i, const unsigned int j, const double aij) {
  if (m == NULL || m->columns == NULL) {
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
  if (m == NULL || m->columns == NULL) {
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

Matrix* matrix_remove_column_vector_at(const Matrix* m, unsigned int column_index) {
  unsigned int i, j;

  if (column_index > m->ncols - 1) {
    errno = EMATRIXINVALIDCOLINDEX;
    return NULL;
  }

  Matrix* n = matrix_create_empty(m->nlins, m->ncols - 1);
  if (n == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  for (i = 0, j = 0; i < m->ncols; i++) {
     if (i == column_index) {
       j = i;
       continue;
     }
     memcpy(n->columns[j].x, m->columns[i].x, m->nlins * sizeof(double));
     j++;
  }

  return n;
}

/**
 * Substitui o vetor coluna de indice i pelo vetor x.
 *
 * A funcao assume a posse de x, ou seja, x nao pode ser liberado pelo
 * chamador. O buffer anterior da coluna e liberado internamente.
 */
void matrix_add_column_vector_at(Matrix* m, unsigned int i, unsigned int xlen, double* x) {
  if (m == NULL || m->columns == NULL) {
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

  Vector vector = {
    .dim = xlen,
    .x = x
  };

  if (m->columns[i].x != NULL) {
    free(m->columns[i].x);
  }

  m->columns[i] = vector;
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

  Matrix* c = matrix_create_empty(a->nlins, a->ncols);
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

  Matrix* c = matrix_create_empty(m->nlins, m->ncols);
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

Matrix* matrix_transpose(const Matrix* m) {
  if (m == NULL) {
    errno = EMATRIXNULL;
    return NULL;
  }

  Vector* line_vectors = matrix_get_line_vectors(m);
  if (line_vectors == NULL) {
    return NULL;
  }

  Matrix* t = matrix_create(m->ncols, m->nlins, line_vectors);

  matrix_free_line_vectors(line_vectors, m->nlins);

  if (t == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  return t;
}

Matrix* matrix_inverse(const Matrix* m) {
  /*
  Entrada: Matrizes L, U, P de tamanho n x n (onde PA = LU)
  Saída: Matriz A_inv (inversa de A)

  1. Inicializar A_inv como uma matriz vazia de tamanho n x n.
  2. Criar a matriz Identidade I de tamanho n x n.

  3. Para cada coluna i de 1 até n:
     a. Obter a i-ésima coluna de I, chamada e_i.
     b. Aplicar a permutação: b = P * e_i (apenas rearranja as linhas do vetor e_i).
     
     c. [Substituição Direta] Resolver L * y = b:
        Para k de 1 até n:
          y[k] = b[k] - (soma de L[k][j] * y[j] para j de 1 até k-1)
          // Nota: L[k][k] é sempre 1 na decomposição LU padrão, 
          // caso contrário, dividiria o resultado por L[k][k].
          
     d. [Substituição Retroativa] Resolver U * x = y:
        Para k de n decrescendo até 1:
          x[k] = (y[k] - (soma de U[k][j] * x[j] para j de k+1 até n)) / U[k][k]
          
     e. Salvar o vetor x como a i-ésima coluna de A_inv (A_inv[:, i] = x).

  4. Retornar A_inv.
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
    // Preservando o tratamento de erro original
    return NULL;
  }

  unsigned int i, j;
  unsigned int n = m->nlins;

  // Inversa a ser retornada
  Matrix* inv = matrix_create_empty(n, n);
  if (inv == NULL) {
    lu_destroy(lu);
    return NULL;
  }

  // Identidade
  Matrix* id = matrix_create_empty(n, n);
  if (id == NULL) {
    matrix_destroy(inv);
    lu_destroy(lu);
    return NULL;
  }

  // Identidade (inv ja vem zerada de matrix_create_empty)
  for (i = 0; i < n; i++) {
    matrix_set(id, i, i, 1.0);
  }

  // Loop principal para cada coluna i da inversa
  for (i = 0; i < n; i++) {
    // e_i = i-ésima coluna de I
    Matrix* ei = matrix_create(n, 1, &id->columns[i]);
    if (ei == NULL) {
      matrix_destroy(inv);
      matrix_destroy(id);
      lu_destroy(lu);
      return NULL;
    }

    // aplica permutação b = P * e_i
    Matrix* b = matrix_mul(lu->p, ei);
    if (b == NULL) {
      matrix_destroy(ei);
      matrix_destroy(inv);
      matrix_destroy(id);
      lu_destroy(lu);
      return NULL;
    }

    // calloc já inicializa com zero, não precisa de loop extra
    double* y = (double*)calloc(n, sizeof(double));
    if (y == NULL) {
      matrix_destroy(b);
      matrix_destroy(ei);
      matrix_destroy(inv);
      matrix_destroy(id);
      lu_destroy(lu);
      errno = EVECTORALLOCATION;
      return NULL;
    }

    // Substituição direta - resolver Ly = b
    for (unsigned int k = 0; k < n; k++) {
      double soma_ly = 0.0;
      for (j = 0; j < k; j++) { // CORRIGIDO: j < k para ir até k-1
        soma_ly += matrix_get(lu->l, k, j) * y[j];
      }
      // CORRIGIDO: O índice da coluna de b é 0, pois b é n x 1
      y[k] = matrix_get(b, k, 0) - soma_ly;
    }

    double* x = (double*)calloc(n, sizeof(double));
    if (x == NULL) {
      free(y);
      matrix_destroy(b);
      matrix_destroy(ei);
      matrix_destroy(inv);
      matrix_destroy(id);
      lu_destroy(lu);
      errno = EVECTORALLOCATION;
      return NULL;
    }

    // Substituição retroativa - resolver Ux = y
    // CORRIGIDO: Cast para int para evitar estouro negativo com unsigned int
    for (int k = (int)n - 1; k >= 0; k--) {
      double soma_ux = 0.0;
      // CORRIGIDO: Varre de k+1 até n, sem embutir a divisão dentro do laço
      for (j = (unsigned int)k + 1; j < n; j++) {
        soma_ux += matrix_get(lu->u, (unsigned int)k, j) * x[j];
      }

      double ukk = matrix_get(lu->u, (unsigned int)k, (unsigned int)k);
      // Proteção contra divisão por zero em matrizes não invertíveis
      if (fabs(ukk) < 1e-9) {
        free(y);
        free(x);
        matrix_destroy(b);
        matrix_destroy(ei);
        matrix_destroy(inv);
        matrix_destroy(id);
        lu_destroy(lu);
        errno = EMATRIXSINGULAR;
        return NULL;
      }
      x[k] = (y[k] - soma_ux) / ukk;
    }

    // Salvar o vetor x como i-ésima coluna da inv.
    // matrix_add_column_vector_at assume a posse de x e libera o buffer
    // anterior da coluna, evitando vazamento.
    matrix_add_column_vector_at(inv, i, n, x);

    // CORRIGIDO: Liberar a memória dos objetos criados NESTA iteração do loop
    free(y);
    matrix_destroy(b);
    matrix_destroy(ei);
  }

  // Limpeza final de memória
  matrix_destroy(id);
  lu_destroy(lu);

  return inv;
}

/**
 * Calcula o sinal da permutacao representada pela matriz P.
 *
 * Retorna 1 para permutacao par, -1 para permutacao impar e 0 em caso de
 * falha de alocacao.
 */
static int matrix_permutation_sign(const Matrix* p) {
  const unsigned int n = p->nlins;

  unsigned int* perm = (unsigned int*)calloc(n, sizeof(unsigned int));
  if (perm == NULL) {
    errno = EMATRIXALLOCATION;
    return 0;
  }

  // P[i][perm[i]] = 1, logo perm[i] e a coluna do unico elemento nao nulo da linha i
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < n; j++) {
      if (matrix_get(p, i, j) != 0.0) {
        perm[i] = j;
        break;
      }
    }
  }

  // O sinal da permutacao e (-1) elevado ao numero de inversoes
  unsigned long inversions = 0;
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = i + 1; j < n; j++) {
      if (perm[i] > perm[j]) {
        inversions++;
      }
    }
  }

  free(perm);

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
  for (unsigned int i = 0; i < n; i++) {
    product *= matrix_get(lu->u, i, i);
  }

  const int sign = matrix_permutation_sign(lu->p);

  lu_destroy(lu);

  if (sign == 0) {
    errno = EMATRIXINVALIDDETERMINANT;
    return NAN;
  }

  return sign * product;
}

void matrix_print(const Matrix* m, const char* matrix_label) {
  if (m == NULL || m->columns == NULL) {
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
