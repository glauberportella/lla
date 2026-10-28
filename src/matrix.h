#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdbool.h>
#include "vector.h"

#define EMATRIXALLOCATION         2000
#define EMATRIXINVALIDADDDIM      2001
#define EMATRIXINVALIDLINEINDEX   2002
#define EMATRIXINVALIDCOLINDEX    2003
#define EMATRIXINVALIDMUL         2004
#define EMATRIXINVALIDDIM         2005
#define EMATRIXSINGULAR           2006
#define EMATRIXINVALIDDETERMINANT 2007
#define EMATRIXNULL               2999

/**
 * A matriz é armazenada em um único bloco contíguo de doubles, no formato
 * column-major: o elemento a(i, j) fica em data[j * nlins + i].
 *
 * Ou seja: a matriz é uma sequência de colunas, mas todas as colunas estão
 * na MESMA alocação (uma matriz n x n é um único array de n*n doubles).
 * Isso mantém a leitura "colunas de vetores" da API sem o custo de uma
 * alocação por coluna e sem dispersar o dado na memória.
 */
typedef struct matrix_rn_t {
  double*       data;  // bloco contíguo column-major, com nlins * ncols doubles
  unsigned int  nlins; // numero de linhas
  unsigned int  ncols; // numero de colunas
} Matrix;

Matrix* matrix_create_empty(unsigned int nlins, unsigned int ncols);
Matrix* matrix_create(unsigned int nlins, unsigned int ncols, Vector columns[]);
Matrix* matrix_create_identity(unsigned int m, unsigned int n);

/**
 * @brief matrix_copy(const Matrix*)
 * Retorna uma copia profunda de m (novo bloco contíguo), ou NULL com errno.
 */
Matrix* matrix_copy(const Matrix* m);

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
 * @brief matrix_get_vector_column_at(const Matrix*, unsigned int)
 * Retorna uma CÓPIA da coluna de indice i de a como um novo Vector.
 *
 * O vetor devolvido deve ser liberado pelo chamador com vector_destroy.
 */
Vector* matrix_get_vector_column_at(const Matrix* a, unsigned int i);

/**
 * @brief matrix_set(Matrix*, const unsigned int, const unsigned int, double)
 * Adiciona um valor na posicao i, j da matriz m
 *
 * @return Matrix* retorna a própria matriz para caso queira concatenar chamadas.
 */
Matrix* matrix_set(Matrix* m, const unsigned int i, const unsigned int j, const double aij);
double matrix_get(const Matrix* m, const unsigned int i, const unsigned int j);

/**
 * @brief matrix_add_column_vector_at(Matrix*, unsigned int, unsigned int, double*)
 * Substitui o conteúdo da coluna de indice i pelos xlen valores de x.
 *
 * ATENÇÃO: diferente das versões anteriores da biblioteca, esta função
 * NÃO assume a posse do buffer x — os valores sao copiados para dentro
 * do bloco da matriz. O chamador continua responsável por liberar x
 * (quando apropriado).
 *
 * Requer xlen igual ao numero de linhas da matriz.
 */
void matrix_add_column_vector_at(Matrix* m, unsigned int i, unsigned int xlen, double* x);
Matrix* matrix_remove_column_vector_at(const Matrix* m, unsigned int column_index);

/**
 * Operações
 *
 * matrix_add, matrix_sub e matrix_mul_scalar devolvem uma nova matriz.
 * matrix_transpose, matrix_inverse e matrix_determinant operam sobre matrizes
 * quadradas, exceto a transposta.
 */
Matrix* matrix_add(const Matrix* a, const Matrix* b);
Matrix* matrix_sub(const Matrix* a, const Matrix* b);
Matrix* matrix_mul(const Matrix* a, const Matrix* b);
Matrix* matrix_mul_scalar(const double scalar, const Matrix* m);
Matrix* matrix_transpose(const Matrix* m);

/**
 * @brief matrix_mul_vector(const Matrix*, const Vector*)
 * Retorna o produto A * x, com A having m linhas e n colunas e x com n
 * coordenadas.
 *
 * @return Vector* o resultado, com m coordenadas, ou NULL com errno definido.
 */
Vector* matrix_mul_vector(const Matrix* m, const Vector* v);

/**
 * @brief matrix_inverse(const Matrix*)
 * Retorna a inversa da matriz quadrada m, resolvendo A * X = I via decomposicao LU.
 *
 * @return Matrix* a inversa, ou NULL com errno = EMATRIXSINGULAR/EMATRIXINVALIDDIM.
 */
Matrix* matrix_inverse(const Matrix* m);

/**
 * @brief matrix_determinant(const Matrix*)
 * Calcula o determinante da matriz quadrada m via decomposicao LU.
 *
 * @return double o determinante, 0.0 para matriz singular, ou NAN com errno
 * definido em caso de erro (matriz nula, nao quadrada ou falha de alocacao).
 */
double matrix_determinant(const Matrix* m);

/**
 * @brief matrix_trace(const Matrix*)
 * Soma da diagonal principal da matriz quadrada m.
 *
 * @return double o traco, ou NAN com errno definido (matriz nula ou nao quadrada).
 */
double matrix_trace(const Matrix* m);

/**
 * @brief matrix_rank(const Matrix*)
 * Calcula o posto da matriz a via decomposicao em valores singulares: o
 * numero de valores singulares acima de eps * max(m, n) * sigma_max.
 *
 * @return unsigned int o posto, ou (unsigned int)-1 com errno definido
 * (matriz nula ou sem linhas/colunas).
 */
unsigned int matrix_rank(const Matrix* m);

/**
 * Saída
 */
void matrix_print(const Matrix* m, const char* matrix_label);
void matrix_println(const Matrix* m, const char* matrix_label);

#endif
