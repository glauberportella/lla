/**
 * Decomposicao LU com pivoteamento parcial
 *
 * A matriz A e decomposta em tres fatores:
 *
 * P*A = L*U
 *
 * onde P e a matriz de permutacao (linhas de A reordenadas), L e triangular
 * inferior com diagonal unitaria e U e triangular superior.
 *
 * Como o pivoteamento parcial reordena as linhas de A, o produto L*U isolado
 * NAO devolve A: e preciso aplicar a permutacao, ou seja, A = P^t * L * U.
 *
 * A = | 4 3 | = LU = | l11 0  | | u11 u12|
 *     | 6 2 |        | l21 l22| | 0   u22|
 *
 * Uma matriz retangular está na sua forma escalonada ou na forma de escada por linhas quando satisfaz as seguintes condições:
 * - Todas as linhas não-nulas estão acima de qualquer linha composta só de zero
 * - O pivô de cada linha está numa coluna à direita do pivô da linha acima.
 * - Todos os elementos de uma coluna abaixo de um pivô são zero.
 *
 *
 * Ideia de algoritmo
 * 1 - Ordenar as linhas referentes aos zeros existentes nas colunas para facilitar achar os pivos
 * 2 - Zerar cada elemento abaixo da diagonal principal: Lij(1) <- Lij - mijL(i+1)j
 */

#ifndef _LU_H_
#define _LU_H_

#include "matrix.h"

typedef struct LU_t {
  Matrix* l; // triangular inferior, diagonal unitaria
  Matrix* u; // triangular superior
  Matrix* p; // permutacao: P*A = L*U
} LU;

LU* lu_create_decomposition(const Matrix* a);
void lu_destroy(LU* lu);

#endif
