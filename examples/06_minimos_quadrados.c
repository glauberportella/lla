/**
 * Exemplo 06 - Minimos quadrados (sistema retangular)
 *
 * Cenario real: um freelancer mede o tempo (horas) que leva para cada
 * tarefa-cliente e quer prever o estimado para uma nova tarefa. Como ha
 * mais medicoes que incognitas, o sistema A*p = y e retangular (m x 2) e
 * a solucao vem dos minimos quadrados via SVD (pseudoinversa).
 *
 * Modelo: tempo = a0 + a1 * tamanho_da_tarefa
 */

#include <stdio.h>
#include <errno.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "lsystem.h"

int main(void) {
  printf("=== 06 - Regressao linear via minimos quadrados ===\n\n");

  // 5 medicoes: tempo (h) x tamanho da tarefa
  double tempo[] = { 0.5, 1.2, 2.1, 2.8, 3.6, 4.1 };      // y
  double tamanho[] = { 1, 2, 4, 5, 7, 8 };               // x0

  // A: coluna 1 = 1 (termo constante), coluna 2 = tamanho
  double c1[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
  double c2[] = { 1.0, 2.0, 4.0, 5.0, 7.0, 8.0 };
  Vector colunas[] = {
    { .dim = 6, .x = c1 },
    { .dim = 6, .x = c2 },
  };
  Matrix* a = matrix_create(6, 2, colunas);
  (void)tamanho;

  Vector* y = vector_create(6, tempo);
  vector_println(y, "y (horas medidos)");

  // lsystem_solve detecta a matriz retangular e usa SVD
  errno = 0;
  Vector* p = lsystem_solve(a, y);
  if (p == NULL) {
    printf("Falha na regressao (errno=%d)\n", errno);
    return 1;
  }

  double a0 = p->x[0];
  double a1 = p->x[1];
  printf("Modelo ajustado: tempo = %.4f + %.4f * tamanho\n", a0, a1);
  printf("Horas fixas por tarefa: %.3f\n", a0);
  printf("Horas por unidade de tamanho: %.3f\n\n", a1);

  // Previsões para duas novas tarefas
  for (int k = 0; k < 2; k++) {
    double tamanho_pred[] = { 3.0, 10.0 };
    printf("Previsao p/ tamanho %.1f: %.2f horas\n",
           tamanho_pred[k], a0 + a1 * tamanho_pred[k]);
  }

  // Residuo: norma do erro do ajuste
  double residuo = 0.0;
  for (unsigned int i = 0; i < 6; i++) {
    double erro = tempo[i] - (a0 + a1 * c2[i]);
    residuo += erro * erro;
  }
  printf("Soma dos erros ao quadrado: %.4f\n", residuo);

  vector_destroy(y);
  vector_destroy(p);
  matrix_destroy(a);

  return 0;
}
