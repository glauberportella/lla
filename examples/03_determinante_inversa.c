/**
 * Exemplo 03 - Determinante, singularidade e inversa
 *
 * Cenario real: engenharia estrutural. Uma estrutura com 3 articulacoes
 * tem equilíbrio estavel se a matriz de rigidez (simplificada aqui) nao
 * for singular. Este exemplo mostra quando o determinante diagnostica a
 * singularidade e como inverter uma matriz invertivel.
 */

#include <stdio.h>
#include <errno.h>
#include <math.h>
#include "matrix.h"
#include "lu.h"

Matrix* construir_stiffness(double k1, double k2, double k3) {
  // Matriz diagonal por blocos simplificada de uma mola com 3 molas:
  // K = | k1  k1  0  |
  //     | k1  k1+k2  k2 |
  //     |  0   k2  k2+k3 |
  double c1[] = { k1, k1, 0.0 };
  double c2[] = { k1, k1 + k2, k2 };
  double c3[] = { 0.0, k2, k2 + k3 };
  Vector colunas[] = {
    { .dim = 3, .x = c1 },
    { .dim = 3, .x = c2 },
    { .dim = 3, .x = c3 },
  };
  return matrix_create(3, 3, colunas);
}

int main(void) {
  printf("=== 03 - Determinante e inversa ===\n\n");

  // Caso 1: estrutura rigida (k1, k2, k3 > 0), matriz invertivel
  Matrix* k_rigida = construir_stiffness(2.0, 3.0, 1.5);
  printf("K rigida:\n");
  matrix_println(k_rigida, "K");

  double det1 = matrix_determinant(k_rigida);
  printf("det(K) = %.4f (diferente de zero: estrutura estavel)\n\n", det1);

  Matrix* k_inv = matrix_inverse(k_rigida);
  matrix_println(k_inv, "K inversa");

  // Verificacao: K * K^-1 deve ser (quase) a identidade
  Matrix* back = matrix_mul(k_rigida, k_inv);
  matrix_println(back, "K * K^-1 (espera ~identidade)");

  // dR/dk2? O traco da inversa e o reciproco da condutividade de mola:
  // fica como comparacao numerica
  printf("Traco da inversa: %.4f\n", matrix_trace(k_inv));
  printf("Traco de K: %.4f\n\n", matrix_trace(k_rigida));

  // Caso 2: mola solta (k3 = 0) - estrutura solta uma direcao, matriz singular
  Matrix* k_solta = construir_stiffness(2.0, 1.0, 0.0);
  double det2 = matrix_determinant(k_solta);
  printf("K solta: det = %.4g (zero: estrutura passa a ser solta)\n", det2);

  // Tentar inverter uma matriz singular volta NULL com errno definido
  errno = 0;
  Matrix* inv_falha = matrix_inverse(k_solta);
  if (inv_falha == NULL) {
    printf("Inversa da matriz singular falhou: errno=%d\n", errno);
    printf("Inversa (espera NULL): %p\n\n", (void*)inv_falha);
  }

  matrix_destroy(k_rigida);
  matrix_destroy(k_inv);
  matrix_destroy(back);
  matrix_destroy(k_solta);
  matrix_destroy(inv_falha);

  return 0;
}
