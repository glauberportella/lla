/**
 * Exemplo 02 - Operacoes com matrizes
 *
 * Cenario real: transformacoes geometricas de um sprite 2D em um jogo.
 * A cada frame o sprite e girado, escalado e transladado. Este exemplo
 * mostra criacao de matrizes, identidade, get/set, insercao e remocao de
 * colunas, adicao, subtracao, produto por escalar, produto de matrizes,
 * transposta, produto por vetor e traco.
 */

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

Matrix* criar_translacao(double tx, double ty) {
  // A translacao em coordenadas homogeneas usa uma matriz 3x3:
  //     | 1 0 tx |
  // T = | 0 1 ty |
  //     | 0 0  1 |
  double c1[] = { 1.0, 0.0, 0.0 };
  double c2[] = { 0.0, 1.0, 0.0 };
  double c3[] = { tx, ty, 1.0 };
  Vector colunas[] = {
    { .dim = 3, .x = c1 },
    { .dim = 3, .x = c2 },
    { .dim = 3, .x = c3 },
  };
  return matrix_create(3, 3, colunas);
}

Matrix* criar_rotacao(double angulo) {
  //     | cos  -sin  0 |
  // R = | sin   cos  0 |
  //     |  0     0   1 |
  double cosseno = cos(angulo);
  double seno = sin(angulo);
  double c1[] = { cosseno, seno, 0.0 };
  double c2[] = { -seno, cosseno, 0.0 };
  double c3[] = { 0.0, 0.0, 1.0 };
  Vector colunas[] = {
    { .dim = 3, .x = c1 },
    { .dim = 3, .x = c2 },
    { .dim = 3, .x = c3 },
  };
  return matrix_create(3, 3, colunas);
}

Matrix* criar_escala(double sx, double sy) {
  //     | sx  0   0 |
  // S = |  0  sy  0 |
  //     |  0   0  1 |
  double c1[] = { sx, 0.0, 0.0 };
  double c2[] = { 0.0, sy, 0.0 };
  double c3[] = { 0.0, 0.0, 1.0 };
  Vector colunas[] = {
    { .dim = 3, .x = c1 },
    { .dim = 3, .x = c2 },
    { .dim = 3, .x = c3 },
  };
  return matrix_create(3, 3, colunas);
}

int main(void) {
  printf("=== 02 - Operacoes com matrizes ===\n\n");

  // Sprite em (4, 2) representado em coordenadas homogeneas
  double pos[] = { 4.0, 2.0, 1.0 };
  Vector sprite[] = { { .dim = 3, .x = pos } };

  Matrix* identidade = matrix_create_identity(3, 3);
  matrix_println(identidade, "Identidade I3");

  Matrix* t = criar_translacao(1.0, 0.5);
  Matrix* r = criar_rotacao(3.14159265358979 / 6.0); // 30 graus
  Matrix* s = criar_escala(2.0, 2.0);

  // Composicao: primeiro escala, depois gira, depois translada
  Matrix* rs = matrix_mul(r, s);
  Matrix* m = matrix_mul(t, rs);

  matrix_println(m, "M = T * R * S");

  // Aplica a transformacao ao sprite
  Vector* nova_pos = matrix_mul_vector(m, sprite);
  vector_println(nova_pos, "Nova posicao do sprite");

  // Traco da matriz (soma da diagonal) e da parte linear
  printf("Traco de M: %.4f\n", matrix_trace(m));
  printf("Traco de R (a soma dos cossenos, espera 2*cos(30o)): %.4f\n\n",
         matrix_trace(r));

  // Produto por escalar e subtracao: m - 2*(m/2) deve virar zero
  Matrix* metade = matrix_mul_scalar(0.5, m);
  Matrix* dobro_da_metade = matrix_mul_scalar(2.0, metade);
  Matrix* resto = matrix_sub(m, dobro_da_metade);
  printf("2*(M/2) e igual a M? (checagem por subtracao)\n");
  matrix_println(resto, "M - 2*(M/2) (espera zeros)");

  // Transposta de M
  Matrix* mt = matrix_transpose(m);
  matrix_println(mt, "M transposta");

  // Validacao de indice e acesso get/set
  if (matrix_validate_index(m, 5, 0) != 0) {
    printf("matrix_get(pos 5,0) e invalido: errno=%d\n", errno);
  }
  printf("M(0,0) antes do set: %.4f\n", matrix_get(m, 0, 0));
  matrix_set(m, 0, 0, 2.0);
  printf("M(0,0) depois do set: %.4f\n\n", matrix_get(m, 0, 0));

  // Adicao de matrizes
  Matrix* soma = matrix_add(m, identidade);
  matrix_println(soma, "M + I");

  // Linhas da matriz como vetores
  Vector* linhas = matrix_get_line_vectors(m);
  vector_println(&linhas[0], "Primeira linha de M");

  // Substituicao de coluna: matrix_add_column_vector_at copia os valores
  // de x para a coluna i e NAO assume a posse do buffer, por isso ele e
  // liberado aqui no fim do programa.
  double* nova_col = (double*)malloc(3 * sizeof(double));
  for (int i = 0; i < 3; i++) {
    nova_col[i] = 7.0;
  }
  printf("Substituindo a coluna 1 por (7,7,7)...\n");
  matrix_add_column_vector_at(m, 1, 3, nova_col);
  matrix_println(m, "M com nova coluna 1");

  printf("Removendo a coluna 2...\n");
  Matrix* m_menos = matrix_remove_column_vector_at(m, 2);
  matrix_println(m_menos, "M 3x2");


  matrix_destroy(identidade);
  matrix_destroy(t);
  matrix_destroy(r);
  matrix_destroy(s);
  matrix_destroy(rs);
  matrix_destroy(m);
  matrix_destroy(metade);
  matrix_destroy(dobro_da_metade);
  matrix_destroy(resto);
  matrix_destroy(mt);
  matrix_destroy(soma);
  matrix_destroy(m_menos);
  vector_destroy(nova_pos);
  free(nova_col);
  for (unsigned int i = 0; i < 3; i++) {
    free(linhas[i].x);
  }
  free(linhas);

  return 0;
}
