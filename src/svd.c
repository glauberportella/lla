#include "svd.h"
#include "matrix.h"
#include "vector.h"
#include <errno.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

int lla_signal(double n) {
  return n >= 0 ? 1 : -1;
}

double* givens_rotation(double f, double g) {
  double c, s;
  
  if (g == 0) {
    c = 1.0;
    s = 0.0;
  } else if (fabs(g) > fabs(f)) {
    double cotan = f / g;
    s = 1 / sqrt(1 + pow(cotan, 2.0));
    c = s * cotan;
  } else {
    double tan = g / f;
    c = 1 / sqrt(1 + pow(tan, 2));
    s = c * tan;
  }
  
  double ret[2] = {c, s};

  return ret;
}

SVDPosProcResult* svd_golub_reinsch(Matrix* a) {
  HouseholderBidiagResult* hbdr = svd_householder_bidiag(a);
  QRBidiagResult* bdr = svd_qr_bidiag_iter(hbdr);
  SVDPosProcResult* ppr = svd_pos_proc(bdr);

  return ppr;
}

HouseholderBidiagResult* svd_householder_bidiag(Matrix* a) {
  unsigned int m = a->nlins;
  unsigned int n = a->ncols;
  
  Matrix* U = matrix_create_identity(m, m);
  if (U == NULL) {
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  Matrix* V = matrix_create_identity(n, n);
  if (V == NULL) {
    matrix_destroy(U);
    errno = EMATRIXALLOCATION;
    return NULL;
  }

  Vector* d = vector_create_empty(n);
  if (d == NULL) {
    matrix_destroy(V);
    matrix_destroy(U);
    errno = EVECTORALLOCATION;
    return NULL;
  }
  
  Vector* e = vector_create_empty(n - 1);
  if (e == NULL) {
    vector_destroy(d);
    matrix_destroy(V);
    matrix_destroy(U);
    errno = EVECTORALLOCATION;
    return NULL;
  }

  unsigned int k;
  
  // vetor de Householder
  Vector* v;

  for (k = 0; k < n - 1; k++) {
    // --- 1. Aniquilar elementos abaixo de A[k,k] (Coluna) ---
    Vector* vector_x = matrix_get_vector_column_at(a, k);
    if (vector_x == NULL) {
      vector_destroy(e);
      vector_destroy(d);
      matrix_destroy(V);
      matrix_destroy(U);
      return NULL;
    }
    double norm_x = vector_norm(vector_x);
    if (norm_x > 0) {
      // PONTO CRÍTICO 3: Evitar cancelamento catastrófico no sinal
      int signal = vector_x->x[0] >= 0 ? 1 : -1;
      double alpha = -signal * norm_x;
      // Construir vetor de Householder v
      v = vector_copy(vector_x);
      v->x[0] = v->x[0] - alpha;
      double norm_v = vector_norm(v);
      if (norm_v > 0) {
        // normaliza vetor refletor
        v = vector_mul(1 / norm_v, v);
        // para fazer vv^T vamos transformar em matriz linha e matriz coluna
        Vector* mcv_columns = {v};
        Matrix* mcv = matrix_create(v->dim, 1, mcv_columns);
        if (mcv == NULL) {
          vector_destroy(vector_x);
          vector_destroy(v);
          vector_destroy(e);
          vector_destroy(d);
          matrix_destroy(V);
          matrix_destroy(U);
          errno = EMATRIXALLOCATION;
          return NULL;
        }
        Matrix* mlv = matrix_transpose(mcv);
        if (mlv == NULL) {
          matrix_destroy(mcv);
          vector_destroy(vector_x);
          vector_destroy(v);
          vector_destroy(e);
          vector_destroy(d);
          matrix_destroy(V);
          matrix_destroy(U);
          errno = EMATRIXALLOCATION;
          return NULL;
        }
        // Aplicar à matriz A à esquerda: A = (I - 2vv^T)A
        Matrix* I = matrix_create_identity(a->nlins, a->ncols);
        if (I == NULL) {
          matrix_destroy(mlv);
          matrix_destroy(mcv);
          vector_destroy(vector_x);
          vector_destroy(v);
          vector_destroy(e);
          vector_destroy(d);
          matrix_destroy(V);
          matrix_destroy(U);
          errno = EMATRIXALLOCATION;
          return NULL;
        }
        Matrix* m2vvT = matrix_mul_scalar(2, matrix_mul(mlv, mcv));
        if (m2vvT == NULL) {
          matrix_destroy(I);
          matrix_destroy(mlv);
          matrix_destroy(mcv);
          vector_destroy(vector_x);
          vector_destroy(v);
          vector_destroy(e);
          vector_destroy(d);
          matrix_destroy(V);
          matrix_destroy(U);
          errno = EMATRIXALLOCATION;
          return NULL;
        }
        a = matrix_mul(matrix_sub(I, m2vvT), a);
        if (a == NULL) {
          matrix_destroy(m2vvT);
          matrix_destroy(I);
          matrix_destroy(mlv);
          matrix_destroy(mcv);
          vector_destroy(vector_x);
          vector_destroy(v);
          vector_destroy(e);
          vector_destroy(d);
          matrix_destroy(V);
          matrix_destroy(U);
          errno = EMATRIXALLOCATION;
          return NULL;
        }
        // Acumular em U: U = U * (I - 2vv^T) à direita nas colunas correspondentes
        U = matrix_mul(U, matrix_sub(I, m2vvT));
        if (U == NULL) {
          matrix_destroy(a);
          matrix_destroy(m2vvT);
          matrix_destroy(I);
          matrix_destroy(mlv);
          matrix_destroy(mcv);
          vector_destroy(vector_x);
          vector_destroy(v);
          vector_destroy(e);
          vector_destroy(d);
          matrix_destroy(V);
          matrix_destroy(U);
          errno = EMATRIXALLOCATION;
          return NULL;
        }

        matrix_destroy(m2vvT);
        matrix_destroy(I);
        matrix_destroy(mlv);
        matrix_destroy(mcv);
        vector_destroy(vector_x);
        vector_destroy(v);
      }
    }

    vector_set(d, k, matrix_get(a, k, k));
    
    // --- 2. Aniquilar elementos à direita de A[k, k+1] (Linha) ---
    if (k < n - 2) {
      Vector* vector_y = matrix_get_vector_column_at(a, k + 1);
      if (vector_y == NULL) {
        vector_destroy(e);
        vector_destroy(d);
        matrix_destroy(V);
        matrix_destroy(U);
        return NULL;
      }

      double norm_y = vector_norm(vector_y);

      if (norm_y > 0) {
        int signal = vector_y->x[0] >= 0 ? 1 : -1;
        double beta = -signal * norm_y;
        v = vector_copy(vector_y);
        v->x[0] = v->x[0] - beta;
        double norm_v = vector_norm(v);
        if (norm_v > 0) {
          v = vector_mul(1 / norm_v, v);
          // para fazer vv^T vamos transformar em matriz linha e matriz coluna
          Vector* mcv_columns = {v};
          Matrix* mcv = matrix_create(v->dim, 1, mcv_columns);
          if (mcv == NULL) {
            vector_destroy(vector_y);
            vector_destroy(v);
            vector_destroy(e);
            vector_destroy(d);
            matrix_destroy(V);
            matrix_destroy(U);
            errno = EMATRIXALLOCATION;
            return NULL;
          }
          Matrix* mlv = matrix_transpose(mcv);
          if (mlv == NULL) {
            matrix_destroy(mcv);
            vector_destroy(vector_y);
            vector_destroy(v);
            vector_destroy(e);
            vector_destroy(d);
            matrix_destroy(V);
            matrix_destroy(U);
            errno = EMATRIXALLOCATION;
            return NULL;
          }
          // Aplicar à matriz A à direita: A = (I - 2vv^T)A
          Matrix* I = matrix_create_identity(a->nlins, a->ncols);
          if (I == NULL) {
            matrix_destroy(mlv);
            matrix_destroy(mcv);
            vector_destroy(vector_y);
            vector_destroy(v);
            vector_destroy(e);
            vector_destroy(d);
            matrix_destroy(V);
            matrix_destroy(U);
            errno = EMATRIXALLOCATION;
            return NULL;
          }
          Matrix* m2vvT = matrix_mul_scalar(2, matrix_mul(mlv, mcv));
          if (m2vvT == NULL) {
            matrix_destroy(I);
            matrix_destroy(mlv);
            matrix_destroy(mcv);
            vector_destroy(vector_y);
            vector_destroy(v);
            vector_destroy(e);
            vector_destroy(d);
            matrix_destroy(V);
            matrix_destroy(U);
            errno = EMATRIXALLOCATION;
            return NULL;
          }
          a = matrix_mul(a, matrix_sub(I, m2vvT));
          if (a == NULL) {
            matrix_destroy(m2vvT);
            matrix_destroy(I);
            matrix_destroy(mlv);
            matrix_destroy(mcv);
            vector_destroy(vector_y);
            vector_destroy(v);
            vector_destroy(e);
            vector_destroy(d);
            matrix_destroy(V);
            matrix_destroy(U);
            errno = EMATRIXALLOCATION;
            return NULL;
          }
          // Acumular em V: V = V * (I - 2vv^T)
          V = matrix_mul(V, matrix_sub(I, m2vvT));
          if (V == NULL) {
            matrix_destroy(m2vvT);
            matrix_destroy(I);
            matrix_destroy(mlv);
            matrix_destroy(mcv);
            vector_destroy(vector_y);
            vector_destroy(v);
            vector_destroy(e);
            vector_destroy(d);
            matrix_destroy(V);
            matrix_destroy(U);
            errno = EMATRIXALLOCATION;
            return NULL;
          }

          matrix_destroy(m2vvT);
          matrix_destroy(I);
          matrix_destroy(mlv);
          matrix_destroy(mcv);
          vector_destroy(vector_y);
          vector_destroy(v);
        }
      }
    }
  }

  HouseholderBidiagResult* result = (HouseholderBidiagResult*)malloc(sizeof(HouseholderBidiagResult));
  if (result == NULL) {
    vector_destroy(e);
    vector_destroy(d);
    matrix_destroy(U);
    matrix_destroy(V);
    return NULL;
  }

  result->e = e;
  result->d = d;
  result->U = U;
  result->V = V;

  return result;
}


QRBidiagResult* svd_qr_bidiag_iter(HouseholderBidiagResult* hbdr) {
  if (hbdr == NULL) {
    errno = ENULL;
    return NULL;
  }

  unsigned int n = hbdr->e->dim;
  double epsilon = 1e-15;
  unsigned int max_iter = 30 * n;
  unsigned int iteration = 0;
  unsigned int p, q;

  while (iteration < max_iter) {
    // PONTO CRÍTICO 1: Critério de Desacoplamento (Splitting)
    // Encontra o menor subbloco de trabalho que vai de p até q.
    q = n;
    while (q > 1 && fabs(hbdr->e->x[q - 2]) < epsilon * (fabs(hbdr->d->x[q - 2]) + fabs(hbdr->d->x[q - 1]))) {
      hbdr->e->x[q - 2] = 0; // convergiu, isola o elemento
      q -= 1;
    }

    // Toda a superdiagonal virou zero. SVD Concluído. Sair do loop
    if (q <= 1) {
      break;
    }

    p = q - 1;
    while (p > 1 && fabs(hbdr->e->x[p -2]) > epsilon * (fabs(hbdr->d->x[p - 2]) + fabs(hbdr->d->x[p - 1]))) {
      p -= 1;
    }

    // Agora trabalhamos estritamente no subbloco de p-1 até q-1
    unsigned int idx_ini = p - 1;
    unsigned int idx_end = q - 1;
    // --- Cálculo do Wilkinson Shift (Deslocamento Implícito) ---
    // Analisa o canto inferior direito 2x2 do bloco atual
    double d1 = hbdr->d->x[idx_end - 1];
    double d2 = hbdr->d->x[idx_end];
    double e1 = hbdr->e->x[idx_end - 1];
    double e0 = (idx_end > 1) ? hbdr->e->x[idx_end - 2] : 0;
    
    // Matriz simétrica B^T * B local
    double f = (d1 * d1) + (e0 * e0);
    double g = d1 * e1;
    double h = (d2 * d2) + (e1 * e1);

    // Deslocamento de Wilkinson puro
    double mu = h - (g * g) / ((f - h) / 2 + lla_signal((f - h) / 2) * sqrt(pow(((f - h) / 2), 2.0) + g * g));
    
    // --- Perseguição do Calombo (Bulge Chasing) ---
    // Configura o primeiro "calombo" fora da bidiagonal
    double x = hbdr->d->x[idx_ini] * hbdr->d->x[idx_ini] - mu;
    double z = hbdr->d->x[idx_ini] * hbdr->e->x[idx_ini];

    for (unsigned k = idx_ini; k < idx_end - 1; k++) {
      // Rotação à direita (Givens) para eliminar z
      double* givens = givens_rotation(x, z);
      double c = givens[0];
      double s = givens[1];
      
      // Aplica rotação nas colunas k e k+1 de d e e
      // Isto gera um calombo indesejado na linha abaixo
      double left_buldge = svd_rotate_columns(hbdr->d, hbdr->e, k, k + 1, c, s);
      
      // Acumula rotação em V
      svd_rotate_matrix_columns(hbdr->V, k, k + 1, c, s);
      
      // Rotação à esquerda (Givens) para eliminar o calombo gerado na linha
      x = hbdr->d->x[k];
      z = left_buldge; // calombo esquerdo
      givens = givens_rotation(x, z);
      c = givens[0];
      s = givens[1];

      double right_buldge = svd_rotate_lines(hbdr->d, hbdr->e, k, k + 1, c, s);

      // Acumula rotação em U
      svd_rotate_matrix_columns(hbdr->U, k, k + 1, c, s);

      // Atualiza coordenadas para o próximo passo da perseguição
      if (k < idx_end - 1) {
        x = hbdr->e->x[k];
        z = right_buldge;
      }    
    }

    iteration++;
  }

  QRBidiagResult* result = (QRBidiagResult*)malloc(sizeof(QRBidiagResult));
  if (result == NULL) {
    errno = ENULL;
    return NULL;
  }
  result->d = hbdr->d;
  result->U = hbdr->U;
  result->V = hbdr->V;

  return result;
}

SVDPosProcResult* svd_pos_proc(QRBidiagResult* bdr) {
  unsigned int n = bdr->d->dim;
  unsigned int i;

  for (i = 0; i < n; i++) {
    if (bdr->d->x[i] < 0) {
      bdr->d->x[i] = -1 * bdr->d->x[i];
      // inverte o vetor correspondente para manter consistência
      for (unsigned int j = 0; j < bdr->V->nlins; j++) {
        bdr->V->columns[i].x[j] = -1 * bdr->V->columns[i].x[j];
      }
    }
  }

  for (i = 1; i < n; i++) {
    double key = bdr->d->x[i];
    unsigned int j = i - 1;
    while (j >= 0 && bdr->d->x[j + 1] < key) {
      double aux = bdr->d->x[j];
      bdr->d->x[j] = bdr->d->x[j + 1];
      bdr->d->x[j + 1] = aux;

      svd_change_columns(bdr->U, j, j + 1);
      svd_change_columns(bdr->V, j, j + 1);
      
      j--;
    }
  }

  Matrix* vt = matrix_transpose(bdr->V);
  if (vt == NULL) {
    return NULL;
  }

  SVDPosProcResult* result = (SVDPosProcResult*)malloc(sizeof(SVDPosProcResult));
  if (result == NULL) {
    return NULL;
  }
  result->U = bdr->U;
  result->sigma = bdr->d;
  result->VT = vt;

  return result;
}

double svd_rotate_columns(Vector* d, Vector* e, unsigned int i, unsigned int j, double c, double s) {
  double d_k = d->x[i];
  double e_k = e->x[i];
  double d_k1 = d->x[j];

  d->x[i] = d_k * c - e_k * s;
  e->x[i] = d_k * s + e_k * c;
  d->x[j] = d_k1 * c;

  double left_buldge = -d_k1 * s;

  return left_buldge;
}

double svd_rotate_lines(Vector* d, Vector* e, unsigned int i, unsigned int j, double c, double s) {
  double e_k = e->x[i];
  double d_k1 = d->x[j];

  e->x[i] = e_k * c - d_k1 * s;
  d->x[j] = e_k * s + d_k1 * c;

  double right_buldge = 0;

  if (i < d->dim - 2) {
    double e_k1 = e->x[j];
    e->x[j] = e_k1 * c;
    right_buldge = -e_k1 * s;
  }

  return right_buldge;
}

void svd_rotate_matrix_columns(Matrix* V, unsigned int i, unsigned int j, double c, double s) {
  unsigned int nlines = V->nlins;
  for (unsigned int l = 0; l < nlines; l++) {
    double val_k = matrix_get(V, l, i);
    double val_k1 = matrix_get(V, l, j);
    // aplica a transformação linear de rotação nas colunas i e j
    matrix_set(V, l, i, val_k * c - val_k1 * s);
    matrix_set(V, l, j, val_k * s + val_k1 * c);
  }
}

void svd_change_columns(Matrix* m, unsigned int i, unsigned int j) {
  for (unsigned int l = 0; l < m->nlins; l++) {
    double aux = matrix_get(m, i, j);
    matrix_set(m, i, j, matrix_get(m, i, j + 1));
    matrix_set(m, i, j + 1, aux);
  }
}
