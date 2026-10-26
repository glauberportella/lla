#include "lsystem.h"
#include "matrix.h"
#include <errno.h>
#include <stddef.h>

Vector* lsystem_solve(const Matrix* a, const Vector* b) {
  /**
   * Resolver Ax = b
   * 
   * Decisão:
   *
   *                       solve(A, b)
   *                             │
   *             ┌───────────────┴───────────────┐
   *       A é Quadrada?                  A é Retangular?
   *             │                               │
   *     ┌───────┴───────┐                       │
   *   [ LU / LUP ]  Determinante = 0?        [ SVD ]
   *  (Super Rápido)       │                  (Garante a resposta
   *                       ▼                  de Mínimos Quadrados)
   *                    [ SVD ]
   *              (Estável, não quebra)
   *
   */

  Vector* sol = {0};

  return sol;
}
