#ifndef _SVD_H_
#define _SVD_H_

#include "matrix.h"
#include "vector.h"

#define ESVDNULL         4000
#define ESVDINVALIDDIM   4001
#define ESVDALLOCATION   4002
#define ESVDNOTCONVERGED 4003

/**
 * Resultado da bidiagonalizacao de Householder.
 *
 * A matriz A e fatorada como A = U * B * V^T, onde B e bidiagonal superior:
 * a diagonal principal fica em d e a superdiagonal em e.
 */
typedef struct householder_bidiag_result_t {
  Vector* d;
  Vector* e;
  Matrix* U;
  Matrix* V;
} HouseholderBidiagResult;

/**
 * Resultado da iteracao QR sobre a matriz bidiagonal.
 *
 * Ao final, U^T * A * V e diagonal e os valores singulares estao em d.
 */
typedef struct qr_bidiag_result_t {
  Vector* d;
  Matrix* U;
  Matrix* V;
} QRBidiagResult;

/**
 * Decomposicao em valores singulares: A = U * diag(sigma) * VT.
 *
 * U tem m linhas e n colunas (m >= n), sigma tem os n valores singulares
 * em ordem decrescente e VT tem n linhas e n colunas.
 */
typedef struct svd_result_t {
  Matrix* U;
  Vector* sigma;
  Matrix* VT;
} SVDResult;

/**
 * @brief svd_golub_reinsch(const Matrix*)
 * Calcula a decomposicao em valores singulares de a pelo algoritmo de
 * Golub-Reinsch (Householder + iteracao QR com deslocamento de Wilkinson).
 *
 * A matriz recebida nao e modificada. Matrizes com mais colunas que linhas
 * sao transpostas internamente.
 *
 * @return SVDResult* a decomposicao, ou NULL com errno definido.
 */
SVDResult* svd_golub_reinsch(const Matrix* a);

/**
 * @brief svd_householder_bidiag(const Matrix*)
 * Reduz a a uma matriz bidiagonal superior B, com A = U * B * V^T.
 *
 * @return HouseholderBidiagResult* o resultado, ou NULL com errno definido.
 */
HouseholderBidiagResult* svd_householder_bidiag(const Matrix* a);

/**
 * @brief svd_qr_bidiag_iter(HouseholderBidiagResult*)
 * Executa a iteracao QR sobre a matriz bidiagonal, acumulando as rotacoes
 * em U e V ate que a superdiagonal se anule.
 *
 * A funcao assume a-posse de d, U e V do resultado recebido: esses ponteiros
 * passam a pertencer ao retorno. O chamador deve liberar o resultado recebido
 * com svd_householder_bidiag_destroy.
 *
 * @return QRBidiagResult* o resultado, ou NULL com errno definido.
 */
QRBidiagResult* svd_qr_bidiag_iter(HouseholderBidiagResult* hbdr);

/**
 * @brief svd_pos_proc(QRBidiagResult*)
 * Normaliza os sinais, ordena os valores singulares em ordem decrescente e
 * monta a matriz VT.
 *
 * A funcao assume a-posse de d e U do resultado recebido. O chamador deve
 * liberar o resultado recebido com svd_qr_bidiag_destroy.
 *
 * @return SVDResult* a decomposicao, ou NULL com errno definido.
 */
SVDResult* svd_pos_proc(QRBidiagResult* bdr);

void svd_destroy(SVDResult* result);
void svd_householder_bidiag_destroy(HouseholderBidiagResult* result);
void svd_qr_bidiag_destroy(QRBidiagResult* result);

/**
 * @brief lla_signal(double)
 * Retorna 1 para valores nao negativos e -1 para valores negativos.
 */
int lla_signal(double v);

/**
 * @brief svd_givens_rotation(double, double, double*, double*)
 * Calcula a rotacao de Givens (c, s) que transforma (f, g) em (r, 0).
 */
void svd_givens_rotation(double f, double g, double* c, double* s);

#endif
