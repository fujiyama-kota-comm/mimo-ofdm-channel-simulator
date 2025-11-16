#include "ofdm.h"

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
 * bit reverse
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
 * 正しい FFT (DIT)
 * X[k] = Σ x[n] e^{-j2π nk/N}
 * ============================================================ */
static void fft_radix2(const complex float *x_in, complex float *X_out, int N) {
  if (!is_power_of_two(N)) {
    fprintf(stderr, "fft_radix2: N=%d is not power of two\n", N);
    exit(1);
  }

  bit_reverse_copy(x_in, X_out, N);

  for (int len = 2; len <= N; len <<= 1) {
    int half = len >> 1;

    for (int start = 0; start < N; start += len) {
      for (int k = 0; k < half; k++) {
        /* ★ここを修正：W = e^{-j 2πk/len} */
        float ang = -2.0f * (float)M_PI * (float)k / (float)len;
        complex float W = cosf(ang) + I * sinf(ang);
        // cos(-θ) + j sin(-θ) = cos θ - j sin θ = e^{-jθ}

        complex float a = X_out[start + k];
        complex float b = X_out[start + k + half] * W;

        X_out[start + k] = a + b;
        X_out[start + k + half] = a - b;
      }
    }
  }
}

/* ============================================================
 * 正しい IFFT (DIT)
 * x[n] = (1/N) Σ X[k] e^{+j2π nk/N}
 * ============================================================ */
static void ifft_radix2(const complex float *X_in, complex float *x_out,
                        int N) {
  if (!is_power_of_two(N)) {
    fprintf(stderr, "ifft_radix2: N=%d is not power of two\n", N);
    exit(1);
  }

  bit_reverse_copy(X_in, x_out, N);

  for (int len = 2; len <= N; len <<= 1) {
    int half = len >> 1;

    for (int start = 0; start < N; start += len) {
      for (int k = 0; k < half; k++) {
        float ang = +2.0f * (float)M_PI * (float)k / (float)len;
        complex float W = cosf(ang) + I * sinf(ang); // e^{+jθ}

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

/* 周波数 → 時間 */
void ofdm_ifft_symbol(const complex float *X_freq, complex float *x_time,
                      const ofdm_cfg_t *cfg) {
  if (!cfg || cfg->nfft <= 0) {
    fprintf(stderr, "ofdm_ifft_symbol: invalid cfg\n");
    exit(1);
  }
  ifft_radix2(X_freq, x_time, cfg->nfft);
}

/* 時間 → 周波数 */
void ofdm_fft_symbol(const complex float *x_time, complex float *X_freq,
                     const ofdm_cfg_t *cfg) {
  if (!cfg || cfg->nfft <= 0) {
    fprintf(stderr, "ofdm_fft_symbol: invalid cfg\n");
    exit(1);
  }
  fft_radix2(x_time, X_freq, cfg->nfft);
}

/* CP付加 */
void ofdm_add_cp(const complex float *x_no_cp, complex float *x_cp,
                 const ofdm_cfg_t *cfg) {
  int N = cfg->nfft;
  int Ncp = cfg->ncp;

  /* 末尾 Ncp を先頭へ */
  for (int i = 0; i < Ncp; i++) {
    x_cp[i] = x_no_cp[N - Ncp + i];
  }
  for (int i = 0; i < N; i++) {
    x_cp[Ncp + i] = x_no_cp[i];
  }
}

/* CP除去 */
void ofdm_remove_cp(const complex float *x_cp, complex float *x_no_cp,
                    const ofdm_cfg_t *cfg) {
  int N = cfg->nfft;
  int Ncp = cfg->ncp;

  for (int i = 0; i < N; i++) {
    x_no_cp[i] = x_cp[Ncp + i];
  }
}

/* OFDM変調 */
void ofdm_modulate_symbol(const complex float *X_freq, complex float *x_out,
                          const ofdm_cfg_t *cfg) {
  int N = cfg->nfft;

  complex float *tmp = malloc(sizeof(complex float) * N);
  if (!tmp) {
    fprintf(stderr, "ofdm_modulate_symbol: malloc failed\n");
    exit(1);
  }

  ofdm_ifft_symbol(X_freq, tmp, cfg);
  ofdm_add_cp(tmp, x_out, cfg);

  free(tmp);
}

/* OFDM逆変調 */
void ofdm_demodulate_symbol(const complex float *x_in, complex float *X_freq,
                            const ofdm_cfg_t *cfg) {
  int N = cfg->nfft;

  complex float *tmp = malloc(sizeof(complex float) * N);
  if (!tmp) {
    fprintf(stderr, "ofdm_demodulate_symbol: malloc failed\n");
    exit(1);
  }

  ofdm_remove_cp(x_in, tmp, cfg);
  ofdm_fft_symbol(tmp, X_freq, cfg);

  free(tmp);
}
