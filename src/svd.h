#ifndef _SVD_H_
#define _SVD_H_

#include "vector.h"
#include "matrix.h"

#define ENULL 4000

typedef struct householder_bidiag_result_t {
  Vector* d;
  Vector* e;
  Matrix* U;
  Matrix* V;
} HouseholderBidiagResult;

typedef struct qr_bidiag_result_t {
  Vector* d;
  Matrix* U;
  Matrix* V;
} QRBidiagResult;

typedef struct svd_posproc_result_t {
  Matrix* U;
  Vector* sigma;
  Matrix* VT;
} SVDPosProcResult;

SVDPosProcResult* svd_golub_reinsch(Matrix* a);
HouseholderBidiagResult* svd_householder_bidiag(Matrix* a);
QRBidiagResult* svd_qr_bidiag_iter(HouseholderBidiagResult* hbdr);
SVDPosProcResult* svd_pos_proc(QRBidiagResult* bdr);

int lla_signal(double v);
double* givens_rotation(double f, double g);
double svd_rotate_columns(Vector* d, Vector* e, unsigned int i, unsigned int j, double c, double s);
double svd_rotate_lines(Vector* d, Vector* e, unsigned int i, unsigned int j, double c, double s);
void svd_rotate_matrix_columns(Matrix* V, unsigned int i, unsigned int j, double c, double s);
void svd_change_columns(Matrix* m, unsigned int i, unsigned int j);

#endif
