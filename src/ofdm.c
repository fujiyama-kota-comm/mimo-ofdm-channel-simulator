/**
 * @file ofdm.c
 * @brief Simple radix-2 FFT / IFFT based OFDM modulator / demodulator.
 *
 *   - FFT:  X[k] = Σ x[n] e^{-j2πnk/N}
 *   - IFFT: x[n] = (1/N) Σ X[k] e^{+j2πnk/N}
 *
 * Cyclic prefix:
 *   - ofdm_add_cp  : appends CP to time-domain symbol
 *   - ofdm_remove_cp: removes CP from received symbol
 */

#include "ofdm.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================
 * Internal helpers
 * ============================================================ */

static int is_power_of_two(int N) { return (N > 0) && ((N & (N - 1)) == 0); }

/* ============================================================
 * Bit-reversed copy
 * ============================================================ */
static void bit_reverse_copy(const complex float *in, complex float *out,
                             int N) {
  int p = 0;
  while ((1 << p) < N)
    p++;

  for (int i = 0; i < N; i++) {
    int r = 0, x = i;
    for (int b = 0; b < p; b++) {
      r = (r << 1) | (x & 1);
      x >>= 1;
    }
    out[r] = in[i];
  }
}

/* ============================================================
 * FFT (DIT, radix-2)
 * X[k] = Σ x[n] e^{-j2πnk/N}
 * ============================================================ */
static void fft_radix2(const complex float *x_in, complex float *X_out, int N) {
  if (!is_power_of_two(N)) {
    fprintf(stderr, "fft_radix2: N=%d is not power of two\n", N);
    exit(EXIT_FAILURE);
  }

  bit_reverse_copy(x_in, X_out, N);

  for (int len = 2; len <= N; len <<= 1) {
    int half = len >> 1;

    for (int start = 0; start < N; start += len) {
      for (int k = 0; k < half; k++) {
        /* Twiddle: W = e^{-j 2πk/len} */
        float ang = -2.0f * (float)M_PI * (float)k / (float)len;
        complex float W = cosf(ang) + I * sinf(ang);

        complex float a = X_out[start + k];
        complex float b = X_out[start + k + half] * W;

        X_out[start + k] = a + b;
        X_out[start + k + half] = a - b;
      }
    }
  }
}

/* ============================================================
 * IFFT (DIT, radix-2)
 * x[n] = (1/N) Σ X[k] e^{+j2πnk/N}
 * ============================================================ */
static void ifft_radix2(const complex float *X_in, complex float *x_out,
                        int N) {
  if (!is_power_of_two(N)) {
    fprintf(stderr, "ifft_radix2: N=%d is not power of two\n", N);
    exit(EXIT_FAILURE);
  }

  bit_reverse_copy(X_in, x_out, N);

  for (int len = 2; len <= N; len <<= 1) {
    int half = len >> 1;

    for (int start = 0; start < N; start += len) {
      for (int k = 0; k < half; k++) {
        float ang = +2.0f * (float)M_PI * (float)k / (float)len;
        complex float W = cosf(ang) + I * sinf(ang); /* e^{+jθ} */

        complex float a = x_out[start + k];
        complex float b = x_out[start + k + half] * W;

        x_out[start + k] = a + b;
        x_out[start + k + half] = a - b;
      }
    }
  }

  float scale = 1.0f / (float)N;
  for (int i = 0; i < N; i++)
    x_out[i] *= scale;
}

/* ============================================================
 * Public OFDM API
 * ============================================================ */

/* Frequency-domain → time-domain (no CP) */
void ofdm_ifft_symbol(const complex float *X_freq, complex float *x_time,
                      const ofdm_cfg_t *cfg) {
  if (!cfg || cfg->nfft <= 0) {
    fprintf(stderr, "ofdm_ifft_symbol: invalid cfg\n");
    exit(EXIT_FAILURE);
  }
  ifft_radix2(X_freq, x_time, cfg->nfft);
}

/* Time-domain → frequency-domain (no CP) */
void ofdm_fft_symbol(const complex float *x_time, complex float *X_freq,
                     const ofdm_cfg_t *cfg) {
  if (!cfg || cfg->nfft <= 0) {
    fprintf(stderr, "ofdm_fft_symbol: invalid cfg\n");
    exit(EXIT_FAILURE);
  }
  fft_radix2(x_time, X_freq, cfg->nfft);
}

/* Add cyclic prefix */
void ofdm_add_cp(const complex float *x_no_cp, complex float *x_cp,
                 const ofdm_cfg_t *cfg) {
  int N = cfg->nfft;
  int Ncp = cfg->ncp;

  /* Copy last Ncp samples of x_no_cp to the head of x_cp */
  for (int i = 0; i < Ncp; i++) {
    x_cp[i] = x_no_cp[N - Ncp + i];
  }
  for (int i = 0; i < N; i++) {
    x_cp[Ncp + i] = x_no_cp[i];
  }
}

/* Remove cyclic prefix */
void ofdm_remove_cp(const complex float *x_cp, complex float *x_no_cp,
                    const ofdm_cfg_t *cfg) {
  int N = cfg->nfft;
  int Ncp = cfg->ncp;

  for (int i = 0; i < N; i++) {
    x_no_cp[i] = x_cp[Ncp + i];
  }
}

/* One OFDM symbol modulation (IFFT + CP) */
void ofdm_modulate_symbol(const complex float *X_freq, complex float *x_out,
                          const ofdm_cfg_t *cfg) {
  int N = cfg->nfft;

  complex float *tmp = malloc(sizeof(complex float) * N);
  if (!tmp) {
    fprintf(stderr, "ofdm_modulate_symbol: malloc failed\n");
    exit(EXIT_FAILURE);
  }

  ofdm_ifft_symbol(X_freq, tmp, cfg);
  ofdm_add_cp(tmp, x_out, cfg);

  free(tmp);
}

/* One OFDM symbol demodulation (CP removal + FFT) */
void ofdm_demodulate_symbol(const complex float *x_in, complex float *X_freq,
                            const ofdm_cfg_t *cfg) {
  int N = cfg->nfft;

  complex float *tmp = malloc(sizeof(complex float) * N);
  if (!tmp) {
    fprintf(stderr, "ofdm_demodulate_symbol: malloc failed\n");
    exit(EXIT_FAILURE);
  }

  ofdm_remove_cp(x_in, tmp, cfg);
  ofdm_fft_symbol(tmp, X_freq, cfg);

  free(tmp);
}
