/**
 * @file ofdm.h
 * @brief Simple OFDM modulation with radix-2 FFT/IFFT.
 *
 * API provided:
 *   - ofdm_ifft_symbol()      : frequency → time (no CP)
 *   - ofdm_fft_symbol()       : time → frequency (no CP)
 *   - ofdm_add_cp()           : append cyclic prefix
 *   - ofdm_remove_cp()        : remove cyclic prefix
 *   - ofdm_modulate_symbol()  : full OFDM modulation (IFFT + CP)
 *   - ofdm_demodulate_symbol(): full OFDM demodulation (CP removal + FFT)
 *
 * FFT size must be a power of two.
 */

#ifndef OFDM_H
#define OFDM_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief OFDM configuration (single symbol) */
typedef struct {
  int nfft; /**< FFT size (must be power of two) */
  int ncp;  /**< Cyclic prefix length (samples)  */
} ofdm_cfg_t;

/* ---------------------- Low-level FFT/IFFT ---------------------- */

/**
 * @brief Compute IFFT of one OFDM symbol (no cyclic prefix).
 *
 * @param X_freq  Input frequency-domain vector (length = nfft)
 * @param x_time  Output time-domain vector   (length = nfft)
 * @param cfg     OFDM configuration (FFT size)
 */
void ofdm_ifft_symbol(const complex float *X_freq, complex float *x_time,
                      const ofdm_cfg_t *cfg);

/**
 * @brief Compute FFT of one OFDM symbol (no cyclic prefix).
 */
void ofdm_fft_symbol(const complex float *x_time, complex float *X_freq,
                     const ofdm_cfg_t *cfg);

/* ---------------------- Cyclic Prefix ---------------------- */

/**
 * @brief Add cyclic prefix.
 *
 * Output length = nfft + ncp
 */
void ofdm_add_cp(const complex float *x_no_cp, complex float *x_cp,
                 const ofdm_cfg_t *cfg);

/**
 * @brief Remove cyclic prefix.
 */
void ofdm_remove_cp(const complex float *x_cp, complex float *x_no_cp,
                    const ofdm_cfg_t *cfg);

/* ---------------------- High-level Mod/Demod ---------------------- */

/**
 * @brief OFDM modulation (IFFT + CP).
 *
 * @param X_freq  Frequency-domain symbols (length = nfft)
 * @param x_out   Time-domain with CP (length = nfft + ncp)
 */
void ofdm_modulate_symbol(const complex float *X_freq, complex float *x_out,
                          const ofdm_cfg_t *cfg);

/**
 * @brief OFDM demodulation (CP removal + FFT).
 */
void ofdm_demodulate_symbol(const complex float *x_in, complex float *X_freq,
                            const ofdm_cfg_t *cfg);

#ifdef __cplusplus
}
#endif

#endif /* OFDM_H */
