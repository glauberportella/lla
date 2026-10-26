/**
 * Exemplo 01 - Operacoes com vetores
 *
 * Cenario real: dois cabos seguram uma placa aposentada em um telhado.
 * Cada cabo exerce uma forca (em Newtons) sobre o ponto de fixacao. Este
 * exemplo mostra como usar a lib para calcular a forca resultante, o
 * trabalho de uma forca ao longo de um deslocamento, normas e o produto
 * vetorial (torque/twist).
 */

#include <stdio.h>
#include <errno.h>
#include <math.h>
#include "vector.h"

int main(void) {
  printf("=== 01 - Operacoes com vetores ===\n\n");

  // Forcas dos cabos em 3D
  double c1_coords[] = { 20.0,  5.0, 50.0 };
  double c2_coords[] = {-10.0, 30.0, 45.0 };
  Vector* cabo1 = vector_create(3, c1_coords);
  Vector* cabo2 = vector_create(3, c2_coords);

  vector_println(cabo1, "F cabo 1 (N)");
  vector_println(cabo2, "F cabo 2 (N)");

  // Forca resultante: soma das forcas
  Vector* resultante = vector_add(cabo1, cabo2);
  vector_println(resultante, "Forca resultante (N)");

  // Modulo das forcas (norma euclidiana)
  printf("Modulo cabo 1: %.2f N\n", vector_norm(cabo1));
  printf("Modulo resultante: %.2f N\n\n", vector_norm(resultante));

  // Deslocamento do ponto de fixacao quando a placa e reposta no lugar
  double d_coords[] = { 0.0, 0.0, -0.15 };
  Vector* deslocamento = vector_create(3, d_coords);

  // Trabalho da resultante: produto interno entre forca e deslocamento
  double trabalho = vector_inner_product(resultante, deslocamento);
  printf("Trabalho da resultante (J): %.4f\n", fabs(trabalho));

  // O cabo 1 nao contribui deslocando perpendicularmente: elevar ao cubo
  // para comparar normas
  printf("Norma p=3 da resultante: %.4f\n", vector_norm_p(resultante, 3));
  printf("Norma infinita da resultante: %.2f (maior componente)\n\n",
         vector_norm_inf(resultante));

  // Braco do cabo em relacao ao eixo: produto vetorial f x d
  double b_coords[] = { 1.0, 0.0, 0.0 };
  Vector braco[] = {
    { .dim = 3, .x = c1_coords },
    { .dim = 3, .x = b_coords },
  };

  Vector* torque = vector_product(3, braco);
  if (torque == NULL) {
    printf("Erro ao calcular produto vetorial (errno=%d)\n", errno);
    return 1;
  }
  vector_println(torque, "Torque F1 x b (N.m)");

  // Escala e subtracao: dobro da resultante menos o cabo 2
  Vector* dobro = vector_mul(2.0, resultante);
  Vector* diferenca = vector_sub(dobro, cabo2);
  vector_println(diferenca, "2*resultante - cabo2");

  // Copia e alteracao de coordenadas
  Vector* ajuste = vector_copy(resultante);
  vector_set(ajuste, 2, ajuste->x[2] - 10.0);
  vector_println(ajuste, "Resultante com componente z ajustada");

  vector_destroy(cabo1);
  vector_destroy(cabo2);
  vector_destroy(resultante);
  vector_destroy(deslocamento);
  vector_destroy(torque);
  vector_destroy(dobro);
  vector_destroy(diferenca);
  vector_destroy(ajuste);

  return 0;
}
