#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdbool.h>
#include "vector.h"

#define EMATRIXALLOCATION       2000
#define EMATRIXINVALIDADDDIM    2001
#define EMATRIXINVALIDLINEINDEX 2002
#define EMATRIXINVALIDCOLINDEX  2003
#define EMATRIXINVALIDMUL       2004
#define EMATRIXINVALIDDIM       2005
#define EMATRIXSINGULAR         2006
#define EMATRIXNULL             2999

typedef struct matrix_rn_t {
  Vector* columns;    // vetores colunas que compoem a matriz
  Vector* lines;      // vetores linhas
  unsigned int nlins; // numero de linhas
  unsigned int ncols; // numero de colunas
} Matrix;

Matrix* matrix_create_empty(unsigned int nlins, unsigned int ncols);
Matrix* matrix_create(unsigned int nlins, unsigned int ncols, Vector columns[]);
void matrix_destroy(Matrix* m);

/**
 * @brief matrix_validate_index(Matrix*, const unsigned int, const unsigned int)
 * Valida se os indices passados estão nos limites do numero de linhas e colunas, respectivamente, da matriz.
 *
 * @return int 
 * Retorna 0 se estão ok, -1 se indice de linha está fora da dimensão e 1 se indice de coluna está fora da dimensão.
 */
int matrix_validate_index(const Matrix* m, const unsigned int i, const unsigned int j);
Vector* matrix_get_line_vectors(const Matrix* a);


/**
 * @brief matrix_set(Matrix*, const unsigned int, const unsigned int, double)
 * Adiciona um valor na posicao i, j da matriz m
 *
 * @return Matrix* retorna a própria matriz para caso queira concatenar chamadas.
 */ 
Matrix* matrix_set(Matrix* m, const unsigned int i, const unsigned int j, const double aij);
double matrix_get(const Matrix* m, const unsigned int i, const unsigned int j);
void matrix_add_column_vector_at(Matrix* m, unsigned int i, unsigned int xlen, double* x);

/**
 * Operações
 */
Matrix* matrix_add(const Matrix* a, const Matrix* b);
Matrix* matrix_sub(const Matrix* a, const Matrix* b);
Matrix* matrix_mul(const Matrix* a, const Matrix* b);
Matrix* matrix_mul_scalar(const double scalar, const Matrix* m);
Matrix* matrix_transpose(const Matrix* m);
Matrix* matrix_inverse(const Matrix* m);
double matrix_determinant(const Matrix* m);

/**
 * Saída
 */
void matrix_print(const Matrix* m, const char* matrix_label);
void matrix_println(const Matrix* m, const char* matrix_label);

#endif
