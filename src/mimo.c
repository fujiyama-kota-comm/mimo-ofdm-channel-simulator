/**
 * @file mimo.c
 * @brief Complex-matrix utilities and linear MIMO equalizers (ZF/MMSE)
 *        for OFDM-based channel simulation.
 *
 * This module implements:
 *   - a minimal complex matrix structure
 *   - matrix operations (transpose, Hermitian, multiply, determinant)
 *   - Gauss–Jordan inversion for small matrices
 *   - per-subcarrier MIMO equalizers:
 *         • Zero-Forcing (ZF)
 *         • Minimum Mean-Square Error (MMSE)
 *
 * Intended for use in MIMO-OFDM channel simulators with
 * 3GPP TDL fading and subcarrier-wise frequency-domain processing.
 */

#include "mimo.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================
 * Internal index helper
 * ============================================================ */
#define IDX(r, c, nc) ((r) * (nc) + (c))

/* ============================================================
 * Matrix allocation / memory management
 * ============================================================ */
mimo_matrix_t *mimo_alloc_matrix(int rows, int cols) {
  mimo_matrix_t *M = malloc(sizeof(mimo_matrix_t));
  if (!M)
    return NULL;
  M->rows = rows;
  M->cols = cols;
  M->data = calloc(rows * cols, sizeof(double complex));
  return M;
}

void mimo_free_matrix(mimo_matrix_t *M) {
  if (!M)
    return;
  free(M->data);
  free(M);
}

void mimo_matrix_copy(const mimo_matrix_t *A, mimo_matrix_t *B) {
  memcpy(B->data, A->data, sizeof(double complex) * A->rows * A->cols);
}

/* ============================================================
 * Basic matrix operations
 * ============================================================ */
void mimo_matrix_set_zero(mimo_matrix_t *A) {
  memset(A->data, 0, sizeof(double complex) * A->rows * A->cols);
}

void mimo_matrix_mul(const mimo_matrix_t *A, const mimo_matrix_t *B,
                     mimo_matrix_t *C) {
  int m = A->rows, n = A->cols, p = B->cols;
  mimo_matrix_set_zero(C);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < p; j++)
      for (int k = 0; k < n; k++)
        C->data[IDX(i, j, p)] += A->data[IDX(i, k, n)] * B->data[IDX(k, j, p)];
}

void mimo_matrix_transpose(const mimo_matrix_t *A, mimo_matrix_t *B) {
  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->cols; j++)
      B->data[IDX(j, i, A->rows)] = A->data[IDX(i, j, A->cols)];
}

void mimo_matrix_conj(const mimo_matrix_t *A, mimo_matrix_t *B) {
  for (int i = 0; i < A->rows * A->cols; i++)
    B->data[i] = conj(A->data[i]);
}

void mimo_matrix_hermite(const mimo_matrix_t *A, mimo_matrix_t *B) {
  /* Hermitian transpose: B = A^H */
  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->cols; j++)
      B->data[IDX(j, i, A->rows)] = conj(A->data[IDX(i, j, A->cols)]);
}

/* ============================================================
 * Determinant:
 *   - Fast paths for 1×1, 2×2, 3×3
 *   - General case: LU decomposition
 * ============================================================ */
double complex mimo_matrix_det(const mimo_matrix_t *A) {
  if (A->rows != A->cols)
    return 0;

  int n = A->rows;

  if (n == 1)
    return A->data[0];

  if (n == 2)
    return A->data[0] * A->data[3] - A->data[1] * A->data[2];

  if (n == 3) {
    double complex a = A->data[0], b = A->data[1], c = A->data[2];
    double complex d = A->data[3], e = A->data[4], f = A->data[5];
    double complex g = A->data[6], h = A->data[7], i = A->data[8];
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
  }

  /* Fallback: LU decomposition */
  mimo_matrix_t *LU = mimo_alloc_matrix(n, n);
  mimo_matrix_copy(A, LU);

  double complex det = 1;

  for (int k = 0; k < n; k++) {
    double complex pivot = LU->data[IDX(k, k, n)];
    if (cabs(pivot) < 1e-12) {
      det = 0;
      break;
    }
    det *= pivot;

    for (int i = k + 1; i < n; i++) {
      LU->data[IDX(i, k, n)] /= pivot;
      for (int j = k + 1; j < n; j++)
        LU->data[IDX(i, j, n)] -=
            LU->data[IDX(i, k, n)] * LU->data[IDX(k, j, n)];
    }
  }

  mimo_free_matrix(LU);
  return det;
}

/* ============================================================
 * Matrix inverse (Gauss–Jordan elimination)
 * ============================================================ */
int mimo_matrix_inverse(const mimo_matrix_t *A, mimo_matrix_t *Ainv) {
  int n = A->rows;
  if (n != A->cols)
    return -1;

  mimo_matrix_t *M = mimo_alloc_matrix(n, n);
  mimo_matrix_copy(A, M);

  /* Identity matrix */
  mimo_matrix_t *Id = mimo_alloc_matrix(n, n);
  for (int i = 0; i < n; i++)
    Id->data[IDX(i, i, n)] = 1.0;

  /* Gauss–Jordan elimination */
  for (int k = 0; k < n; k++) {

    double complex pivot = M->data[IDX(k, k, n)];
    if (cabs(pivot) < 1e-12) {
      mimo_free_matrix(M);
      mimo_free_matrix(Id);
      return -2; /* singular matrix */
    }

    double complex inv_pivot = 1.0 / pivot;

    /* Normalize pivot row */
    for (int j = 0; j < n; j++) {
      M->data[IDX(k, j, n)] *= inv_pivot;
      Id->data[IDX(k, j, n)] *= inv_pivot;
    }

    /* Eliminate column k */
    for (int i = 0; i < n; i++) {
      if (i == k)
        continue;

      double complex factor = M->data[IDX(i, k, n)];

      for (int j = 0; j < n; j++) {
        M->data[IDX(i, j, n)] -= factor * M->data[IDX(k, j, n)];
        Id->data[IDX(i, j, n)] -= factor * Id->data[IDX(k, j, n)];
      }
    }
  }

  /* Output inverse */
  mimo_matrix_copy(Id, Ainv);

  mimo_free_matrix(M);
  mimo_free_matrix(Id);
  return 0;
}

/* ============================================================
 * Zero-Forcing equalizer: X = (H^H H)^(-1) H^H Y
 * ============================================================ */
void mimo_equalize_zf(mimo_matrix_t **H_f, mimo_matrix_t **Y_f,
                      mimo_matrix_t **X_hat, int Nsub, int Nrx, int Ntx) {
  for (int sc = 0; sc < Nsub; sc++) {

    mimo_matrix_t *H = H_f[sc];
    mimo_matrix_t *Y = Y_f[sc];
    mimo_matrix_t *Xo = X_hat[sc];

    mimo_matrix_t *Hh = mimo_alloc_matrix(Ntx, Nrx);
    mimo_matrix_t *HHH = mimo_alloc_matrix(Ntx, Ntx);
    mimo_matrix_t *HHH_inv = mimo_alloc_matrix(Ntx, Ntx);
    mimo_matrix_t *W = mimo_alloc_matrix(Ntx, Nrx);

    mimo_matrix_hermite(H, Hh);
    mimo_matrix_mul(Hh, H, HHH);
    mimo_matrix_inverse(HHH, HHH_inv);
    mimo_matrix_mul(HHH_inv, Hh, W);

    mimo_matrix_mul(W, Y, Xo);

    mimo_free_matrix(Hh);
    mimo_free_matrix(HHH);
    mimo_free_matrix(HHH_inv);
    mimo_free_matrix(W);
  }
}

/* ============================================================
 * MMSE equalizer:
 *       X = (H^H H + σ² I)^(-1) H^H Y
 * ============================================================ */
void mimo_equalize_mmse(mimo_matrix_t **H_f, mimo_matrix_t **Y_f,
                        mimo_matrix_t **X_hat, int Nsub, int Nrx, int Ntx,
                        double noise_var, double tx_power) {
  for (int sc = 0; sc < Nsub; sc++) {

    mimo_matrix_t *H = H_f[sc];
    mimo_matrix_t *Y = Y_f[sc];
    mimo_matrix_t *Xo = X_hat[sc];

    mimo_matrix_t *Hh = mimo_alloc_matrix(Ntx, Nrx);
    mimo_matrix_t *HHH = mimo_alloc_matrix(Ntx, Ntx);
    mimo_matrix_t *R = mimo_alloc_matrix(Ntx, Ntx);
    mimo_matrix_t *Rinv = mimo_alloc_matrix(Ntx, Ntx);
    mimo_matrix_t *W = mimo_alloc_matrix(Ntx, Nrx);

    mimo_matrix_hermite(H, Hh);
    mimo_matrix_mul(Hh, H, HHH);

    /* Add σ² I to H^H H */
    mimo_matrix_copy(HHH, R);
    for (int i = 0; i < Ntx; i++)
      R->data[IDX(i, i, Ntx)] += noise_var;

    mimo_matrix_inverse(R, Rinv);
    mimo_matrix_mul(Rinv, Hh, W);

    mimo_matrix_mul(W, Y, Xo);

    mimo_free_matrix(Hh);
    mimo_free_matrix(HHH);
    mimo_free_matrix(R);
    mimo_free_matrix(Rinv);
    mimo_free_matrix(W);
  }
}
