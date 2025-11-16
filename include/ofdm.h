#ifndef OFDM_H
#define OFDM_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * OFDM configuration
 * ============================================================ */
typedef struct {
  int nfft; /* FFT size (must be power of 2) */
  int ncp;  /* Cyclic prefix length (samples) */
} ofdm_cfg_t;

/* ------------------------------------------------------------
 * Low-level: IFFT / FFT (1 symbol, no CP)
 *   X_freq: length = cfg->nfft
 *   x_time: length = cfg->nfft
 * ---------------------------------------------------------- */
void ofdm_ifft_symbol(const complex float *X_freq, complex float *x_time,
                      const ofdm_cfg_t *cfg);

void ofdm_fft_symbol(const complex float *x_time, complex float *X_freq,
                     const ofdm_cfg_t *cfg);

/* ------------------------------------------------------------
 * CP add / remove (1 symbol)
 * ---------------------------------------------------------- */
void ofdm_add_cp(const complex float *x_no_cp, complex float *x_cp,
                 const ofdm_cfg_t *cfg);

void ofdm_remove_cp(const complex float *x_cp, complex float *x_no_cp,
                    const ofdm_cfg_t *cfg);

/* ------------------------------------------------------------
 * High-level OFDM modulation/demodulation
 * ---------------------------------------------------------- */
void ofdm_modulate_symbol(const complex float *X_freq, complex float *x_out,
                          const ofdm_cfg_t *cfg);

void ofdm_demodulate_symbol(const complex float *x_in, complex float *X_freq,
                            const ofdm_cfg_t *cfg);

#ifdef __cplusplus
}
#endif

#endif /* OFDM_H */
