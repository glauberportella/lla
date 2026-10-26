#ifndef _LSYSTEM_H_
#define _LSYSTEM_H_

#include "matrix.h"
#include "vector.h"

#define ELSYSTEMNULL        3000
#define ELSYSTEMINVALIDDIM  3001

/**
 * @brief lsystem_solve(const Matrix*, const Vector*)
 * Resolve o sistema linear A * x = b.
 *
 * Para matrizes quadradas o sistema e resolvido via decomposicao LU com
 * pivoteamento parcial; se a matriz e singular, ou o sistema e retangular
 * (excesso/pouco definido), a solucao no sentido dos minimos quadrados e
 * obtida via decomposicao em valores singulares (pseudoinversa), com
 * norma minima entre as solucoes.
 *
 * @return Vector* a solucao x de dimensao a->ncols, ou NULL com errno definido.
 */
Vector* lsystem_solve(const Matrix* a, const Vector* b);

#endif
