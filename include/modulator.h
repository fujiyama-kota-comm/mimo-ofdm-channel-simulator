#ifndef MODULATOR_H
#define MODULATOR_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------
 * Modulation Types
 * ----------------------------- */
typedef enum {
  MOD_BPSK = 1,
  MOD_QPSK = 2,
  MOD_16QAM = 4,
  MOD_64QAM = 6,
  MOD_256QAM = 8
} modulation_t;

/* -----------------------------
 *  Mapping / Demapping
 * ----------------------------- */
void mod_map(const int *bits, complex float *symbols, int n_symbols,
             modulation_t mod);

void mod_demod(const complex float *symbols, int *bits, int n_symbols,
               modulation_t mod);

#ifdef __cplusplus
}
#endif

#endif /* MODULATOR_H */
