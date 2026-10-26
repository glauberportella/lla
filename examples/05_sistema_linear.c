/**
 * Exemplo 05 - Sistema linear quadrado (lsystem_solve)
 *
 * Cenario real: planejamento de dieta. Um nutricionista monta uma refeicao
 * com 3 alimentos para acertar exatamente as quantidades de proteina,
 * carboidrato e gordura. Os tres balancos formam um sistema 3x3:
 *
 *      [p1]        [gramas de alimento 1]
 *  A * [p2]    =   [gramas de alimento 2]
 *      [p3]        [gramas de alimento 3]
 */

#include <stdio.h>
#include <errno.h>
#include "matrix.h"
#include "vector.h"
#include "lsystem.h"

int main(void) {
  printf("=== 05 - Sistema linear 3x3 (dieta) ===\n\n");

  // Nutrientes por 100g de cada alimento (colunas = alimentos)
  //      frango  arroz  brócolis
  // prot    31     2.5      2.8
  // carb     0    28       5.0
  // gord     3.5   0.4      0.4
  double prot[] = { 31.0, 2.5, 2.8 };
  double carb[] = { 0.0, 28.0, 5.0 };
  double gord[] = { 3.5, 0.4, 0.4 };
  Vector colunas[] = {
    { .dim = 3, .x = prot },
    { .dim = 3, .x = carb },
    { .dim = 3, .x = gord },
  };
  Matrix* a = matrix_create(3, 3, colunas);
  matrix_println(a, "A (nutrientes por coluna)");

  // Meta: ~70 de proteina, 61 de carboidrato e 8 de gordura
  // (em unidades da tabela, por 100g)
  double meta[] = { 65.5, 61.4, 16.0 };
  Vector* b = vector_create(3, meta);
  vector_println(b, "b (metas de macro em unidades)");

  // Solucao: lsystem_solve usa LU com pivoteamento parcial
  errno = 0;
  Vector* x = lsystem_solve(a, b);
  if (x == NULL) {
    printf("Falha ao resolver o sistema (errno=%d)\n", errno);
    return 1;
  }

  printf("Quantidades (x esta em porcoes de 100g)\n");
  vector_println(x, "x");
  printf("Frango:   %.2f x 100g\n", x->x[0]);
  printf("Arroz:    %.2f x 100g\n", x->x[1]);
  printf("Brocolis: %.2f x 100g\n\n", x->x[2]);

  // Verificacao: A * x deve reproduzir as metas b
  Vector* verificao = matrix_mul_vector(a, x);
  vector_println(verificao, "A*x (espera b)");

  vector_destroy(b);
  vector_destroy(x);
  vector_destroy(verificao);
  matrix_destroy(a);

  return 0;
}
