/**
 * @file modulator.h
 * @brief Gray-mapped digital modulations (BPSK / QPSK / 16/64/256QAM).
 *
 * Bit order:
 *   - BPSK: 1 bit per symbol
 *   - QPSK: b0 = I(MSB), b1 = Q(LSB)
 *   - 16QAM: 2 bits for I, then 2 bits for Q (MSB first)
 *   - 64QAM: 3 bits for I, then 3 bits for Q
 *   - 256QAM: 4 bits for I, then 4 bits for Q
 *
 * Average symbol energy Es is normalized to 1.
 */

#ifndef MODULATOR_H
#define MODULATOR_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Supported modulation orders */
typedef enum {
  MOD_BPSK = 1,  /**< 1 bit/sym  */
  MOD_QPSK = 2,  /**< 2 bits/sym */
  MOD_16QAM = 4, /**< 4 bits/sym */
  MOD_64QAM = 6, /**< 6 bits/sym */
  MOD_256QAM = 8 /**< 8 bits/sym */
} modulation_t;

/**
 * @brief Map binary bits → complex symbols.
 *
 * @param bits        Input bits (length = n_symbols * bits_per_symbol)
 * @param symbols     Output complex symbols (length = n_symbols)
 * @param n_symbols   Number of symbols
 * @param mod         Modulation scheme
 */
void mod_map(const int *bits, complex float *symbols, int n_symbols,
             modulation_t mod);

/**
 * @brief Hard-decision demapper: symbols → bits.
 *
 * @param symbols     Input symbols (equalized)
 * @param bits        Output bits
 * @param n_symbols   Number of symbols
 * @param mod         Modulation scheme
 */
void mod_demod(const complex float *symbols, int *bits, int n_symbols,
               modulation_t mod);

#ifdef __cplusplus
}
#endif

#endif /* MODULATOR_H */
