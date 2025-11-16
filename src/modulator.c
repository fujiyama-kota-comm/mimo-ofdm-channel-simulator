#include "modulator.h"
#include <math.h>
#include <stdio.h>

/* ============================================================
 * Helper: clamp
 * ============================================================ */
static inline int clamp_int(int x, int min, int max) {
  return x < min ? min : (x > max ? max : x);
}

/* ============================================================
 * Gray encode / decode helpers
 * ============================================================ */
static inline int gray_encode(int x) { return x ^ (x >> 1); }

static int gray_decode(int g) {
  int x = g;
  while (g >>= 1) {
    x ^= g;
  }
  return x;
}

/* ============================================================
 * BPSK / QPSK mapping
 * ============================================================ */
static complex float map_bpsk(int b) {
  /* Bits: 0 -> +1, 1 -> -1 (Es = 1) */
  return (b ? -1.0f : 1.0f) + 0.0f * I;
}

static complex float map_qpsk(int b0, int b1) {
  /* Gray: b0 = MSB (I), b1 = LSB (Q)
     00: (+,+), 01:(+,-), 11:(-,-), 10:(-,+)
     Es = 1 with 1/sqrt(2) scaling
  */
  float inv_sqrt2 = 1.0f / sqrtf(2.0f);
  float x = (b0 ? -1.0f : 1.0f) * inv_sqrt2;
  float y = (b1 ? -1.0f : 1.0f) * inv_sqrt2;
  return x + y * I;
}

/* ============================================================
 * PAM Gray mapping for one axis (16/64/256QAM)
 * ============================================================ */
/* Map k bits (MSB first) to normalized PAM level with Gray coding.
 * L = 2^k levels, values: ±1, ±3, ..., ±(L-1)
 * Normalization such that E{|S|^2} = 1 over 2D QAM.
 */
static float pam_gray_map_axis(const int *bits, int k) {
  int L = 1 << k;

  /* Build binary index from bits: bits[0] is MSB */
  int b = 0;
  for (int i = 0; i < k; i++) {
    b = (b << 1) | (bits[i] & 1);
  }

  /* Binary -> Gray index on axis */
  int g = gray_encode(b);

  /* Gray index g in [0, L-1] -> PAM level ±1,±3,... */
  int pam = 2 * g - (L - 1); /* e.g., L=4: g=0..3 -> -3,-1,1,3 */

  /* Normalization: Es_total = 1
     E[pam^2] = (L^2 - 1)/3
     For QAM, Es = 2 * E[(pam/norm)^2] = 1
     -> norm^2 = (2/3)*(L^2 - 1)
  */
  float norm = sqrtf((2.0f / 3.0f) * (float)(L * L - 1));

  return pam / norm;
}

/* Inverse: demap normalized value on one axis to k Gray bits */
static void pam_gray_demap_axis(float val, int *bits, int k) {
  int L = 1 << k;
  float norm = sqrtf((2.0f / 3.0f) * (float)(L * L - 1));

  /* Approximate original PAM level: pam_est ≈ pam ∈ {...,-3,-1,1,3,...} */
  float pam_est = val * norm;

  /* Index on axis: g_hat = (pam + (L-1))/2 in [0, L-1] */
  float tmp = (pam_est + (float)(L - 1)) / 2.0f;
  int g_hat = clamp_int((int)lroundf(tmp), 0, L - 1);

  /* Gray -> binary index */
  int b = gray_decode(g_hat);

  /* binary index -> bits[0..k-1] (MSB first) */
  for (int i = 0; i < k; i++) {
    int shift = k - 1 - i;
    bits[i] = (b >> shift) & 1;
  }
}

/* ============================================================
 * Mapping main function
 * ============================================================ */
void mod_map(const int *bits, complex float *symbols, int n_symbols,
             modulation_t mod) {
  int idx = 0;

  for (int n = 0; n < n_symbols; n++) {

    switch (mod) {

    case MOD_BPSK: {
      int b = bits[idx++];
      symbols[n] = map_bpsk(b);
      break;
    }

    case MOD_QPSK: {
      int b0 = bits[idx++];
      int b1 = bits[idx++];
      symbols[n] = map_qpsk(b0, b1);
      break;
    }

    case MOD_16QAM: {
      /* 4 bits per symbol: 2 for I, 2 for Q */
      float Ire = pam_gray_map_axis(&bits[idx], 2);
      idx += 2;
      float Qim = pam_gray_map_axis(&bits[idx], 2);
      idx += 2;
      symbols[n] = Ire + Qim * I;
      break;
    }

    case MOD_64QAM: {
      /* 6 bits per symbol: 3 for I, 3 for Q */
      float Ire = pam_gray_map_axis(&bits[idx], 3);
      idx += 3;
      float Qim = pam_gray_map_axis(&bits[idx], 3);
      idx += 3;
      symbols[n] = Ire + Qim * I;
      break;
    }

    case MOD_256QAM: {
      /* 8 bits per symbol: 4 for I, 4 for Q */
      float Ire = pam_gray_map_axis(&bits[idx], 4);
      idx += 4;
      float Qim = pam_gray_map_axis(&bits[idx], 4);
      idx += 4;
      symbols[n] = Ire + Qim * I;
      break;
    }

    default:
      fprintf(stderr, "mod_map: Unsupported modulation: %d\n", mod);
      return;
    }
  }
}

/* ============================================================
 * Inverse mapping (hard decision)
 * ============================================================ */
void mod_demod(const complex float *symbols, int *bits, int n_symbols,
               modulation_t mod) {
  int idx = 0;

  for (int n = 0; n < n_symbols; n++) {

    float Ire = crealf(symbols[n]);
    float Qim = cimagf(symbols[n]);

    switch (mod) {

    case MOD_BPSK:
      bits[idx++] = (Ire < 0.0f) ? 1 : 0;
      break;

    case MOD_QPSK:
      bits[idx++] = (Ire < 0.0f) ? 1 : 0;
      bits[idx++] = (Qim < 0.0f) ? 1 : 0;
      break;

    case MOD_16QAM: {
      /* 2 bits for I, 2 bits for Q */
      pam_gray_demap_axis(Ire, &bits[idx], 2);
      idx += 2;
      pam_gray_demap_axis(Qim, &bits[idx], 2);
      idx += 2;
      break;
    }

    case MOD_64QAM: {
      pam_gray_demap_axis(Ire, &bits[idx], 3);
      idx += 3;
      pam_gray_demap_axis(Qim, &bits[idx], 3);
      idx += 3;
      break;
    }

    case MOD_256QAM: {
      pam_gray_demap_axis(Ire, &bits[idx], 4);
      idx += 4;
      pam_gray_demap_axis(Qim, &bits[idx], 4);
      idx += 4;
      break;
    }

    default:
      fprintf(stderr, "mod_demod: Unsupported modulation: %d\n", mod);
      return;
    }
  }
}
