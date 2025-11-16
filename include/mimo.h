/**
 * @file mimo.h
 * @brief Linear-algebra utilities and MIMO equalizers for OFDM channel
 * simulation.
 *
 * This module provides:
 *   - a lightweight complex-matrix structure (`mimo_matrix_t`)
 *   - basic matrix operations (transpose, Hermitian, multiplication,
 * determinant)
 *   - Gauss–Jordan matrix inversion
 *   - MIMO linear equalizers (ZF / MMSE) applied per OFDM subcarrier
 *
 * Designed for use in MIMO-OFDM channel simulators with
 * 3GPP TDL fading and subcarrier-wise frequency-domain processing.
 */

#ifndef MIMO_H
#define MIMO_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * Matrix Structure
 * ============================================================ */
typedef struct {
  int rows;
  int cols;
  double complex *data; /* length = rows * cols */
} mimo_matrix_t;

/* ============================================================
 * Memory management
 * ============================================================ */
mimo_matrix_t *mimo_alloc_matrix(int rows, int cols);
void mimo_free_matrix(mimo_matrix_t *M);
void mimo_matrix_copy(const mimo_matrix_t *A, mimo_matrix_t *B);

/* ============================================================
 * Basic matrix operations
 * ============================================================ */
void mimo_matrix_set_zero(mimo_matrix_t *A);
void mimo_matrix_mul(const mimo_matrix_t *A, const mimo_matrix_t *B,
                     mimo_matrix_t *C);

void mimo_matrix_transpose(const mimo_matrix_t *A, mimo_matrix_t *B);
void mimo_matrix_conj(const mimo_matrix_t *A, mimo_matrix_t *B);
void mimo_matrix_hermite(const mimo_matrix_t *A, mimo_matrix_t *B);

/* Matrix inverse via Gauss–Jordan.
 * Returns 0 on success, negative on singular/invalid matrix.
 */
int mimo_matrix_inverse(const mimo_matrix_t *A, mimo_matrix_t *Ainv);

/* Determinant (optimized for 1×1, 2×2, 3×3; fallback to LU) */
double complex mimo_matrix_det(const mimo_matrix_t *A);

/* ============================================================
 * MIMO Equalizers (per OFDM subcarrier)
 * ============================================================ */

/* Zero-Forcing equalizer:
 *   X_hat = (H^H H)^(-1) H^H Y
 */
void mimo_equalize_zf(mimo_matrix_t **H_f,   /* Nsub × (Nrx × Ntx) */
                      mimo_matrix_t **Y_f,   /* Nsub × (Nrx × 1)    */
                      mimo_matrix_t **X_hat, /* Nsub × (Ntx × 1)    */
                      int Nsub, int Nrx, int Ntx);

/* MMSE equalizer:
 *   X_hat = (H^H H + σ² I)^(-1) H^H Y
 *
 * noise_var : noise variance per complex dimension
 * tx_power  : average Es per spatial stream
 */
void mimo_equalize_mmse(mimo_matrix_t **H_f, mimo_matrix_t **Y_f,
                        mimo_matrix_t **X_hat, int Nsub, int Nrx, int Ntx,
                        double noise_var, double tx_power);

#ifdef __cplusplus
}
#endif

#endif /* MIMO_H */
