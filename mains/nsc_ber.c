#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN32
#include <direct.h> /* Windows: _mkdir */
#else
#include <sys/stat.h> /* POSIX mkdir */
#include <sys/types.h>
#endif

#include "nsc_decoder.h"
#include "nsc_encoder.h"

#define PI 3.14159265358979323846

/* ============================================================================
 *  Monte-Carlo BER Simulation for NSC Code (Rate 1/2)
 * ============================================================================
 *
 *  This program evaluates the BER performance of the non-systematic
 *  convolutional (NSC) code using:
 *
 *      • AWGN channel (BPSK modulation)
 *      • Soft-decision Viterbi decoding (LLR)
 *      • Hard-decision Viterbi decoding (Hamming metric)
 *      • Comparison with uncoded BPSK theoretical BER
 *
 *  Output CSV:
 *      results/nsc_ber_data.csv
 *
 *  Format:
 *      EbN0_dB, BER_soft, BER_hard, BER_bpsk
 *
 *  Typical usage:
 *      $ make
 *      $ ./nsc_ber
 *      $ python python/plot_nsc_ber.py
 *
 * ============================================================================
 */

/* ---------------------------------------------------------------------------
 * Simulation configuration
 * -------------------------------------------------------------------------- */
int trials = 100000; /* Number of Monte-Carlo frames */
double EbN0_min = 0.0;
double EbN0_max = 10.0;
double EbN0_step = 1.0;

/* ============================================================================
 *  Gaussian noise generator: Box–Muller
 * ============================================================================
 *  Generates N(0,1) random variables for AWGN simulation.
 * ========================================================================== */
double randn() {
  double u1 = (rand() + 1.0) / (RAND_MAX + 2.0);
  double u2 = (rand() + 1.0) / (RAND_MAX + 2.0);
  return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

/* ============================================================================
 *  Theoretical BER of Uncoded BPSK over AWGN
 * ============================================================================
 *      BER = Q( sqrt(2 * Eb/N0) )
 *          = 0.5 * erfc( sqrt(Eb/N0) )
 *
 *  Used for comparison with coded performance.
 * ========================================================================== */
double bpsk_ber(double EbN0_linear) { return 0.5 * erfc(sqrt(EbN0_linear)); }

/* ============================================================================
 *  MAIN: BER Simulation for NSC (Viterbi)
 * ============================================================================
 *  - NSC encoder: rate-1/2, terminated (tail bits = {0,0})
 *  - Soft Viterbi and hard Viterbi decoding
 *  - AWGN channel with BPSK mapping
 *
 *  Results stored in results/nsc_ber_data.csv
 * ========================================================================== */
int main() {

  /* ----------------------------------------------------------------------
   * Create results/ directory (ignore if already exists)
   * ------------------------------------------------------------------ */
#ifdef _WIN32
  _mkdir("results");
#else
  mkdir("results", 0777);
#endif

  /* ----------------------------------------------------------------------
   * NSC code configuration
   * ----------------------------------------------------------------------
   *   K = nsc_info_len
   *   N = 2*(K + tail_len)
   *
   *   tail_len = 2 → terminated code
   * ------------------------------------------------------------------ */
  nsc_info_len = 100;           /* K: information bits */
  nsc_code_len = 2 * (100 + 2); /* N: encoded length   */

  int K = nsc_info_len;
  int N = nsc_code_len;

  /* ----------------------------------------------------------------------
   * Open CSV output
   * ------------------------------------------------------------------ */
  FILE *fp = fopen("results/nsc_ber_data.csv", "w");
  if (!fp) {
    fprintf(stderr, "Cannot open results/nsc_ber_data.csv\n");
    return -1;
  }
  fprintf(fp, "EbN0_dB,BER_soft,BER_hard,BER_bpsk\n");

  /* ----------------------------------------------------------------------
   * Allocate buffers
   * ------------------------------------------------------------------ */
  int *data = malloc(sizeof(int) * K);
  int *code = malloc(sizeof(int) * N);
  double *LLR = malloc(sizeof(double) * N);
  int *rx_bits = malloc(sizeof(int) * N);
  int *info_soft = malloc(sizeof(int) * K);
  int *info_hard = malloc(sizeof(int) * K);
  int *code_hat = malloc(sizeof(int) * N); /* optional re-encoded bits */

  if (!data || !code || !LLR || !rx_bits || !info_soft || !info_hard ||
      !code_hat) {
    fprintf(stderr, "Memory allocation failed\n");
    return -1;
  }

  /* RNG seed */
  srand((unsigned)time(NULL));

  printf("EbN0_dB, BER_soft, BER_hard, BER_bpsk\n");

  /* ==========================================================================
   *  Eb/N0 loop
   * ======================================================================== */
  for (double EbN0_dB = EbN0_min; EbN0_dB <= EbN0_max; EbN0_dB += EbN0_step) {

    double EbN0 = pow(10.0, EbN0_dB / 10.0);
    double R = 0.5;                      /* Rate of NSC code */
    double var = 1.0 / (2.0 * R * EbN0); /* Noise variance σ² */
    double sigma = sqrt(var);

    long total_bits = 0;
    long error_soft = 0;
    long error_hard = 0;

    /* ------------------------------------------------------------------
     * Monte-Carlo loop
     * ------------------------------------------------------------------ */
    for (int t = 0; t < trials; t++) {

      /* --- 1) Generate random information bits --- */
      for (int i = 0; i < K; i++)
        data[i] = rand() & 1;

      /* --- 2) Encode (rate-1/2 NSC, terminated) --- */
      nsc_encode_r05(data, code);

      /* --- 3) AWGN channel with BPSK mapping --- */
      for (int i = 0; i < N; i++) {

        double s = (code[i] == 0 ? +1.0 : -1.0); /* BPSK */

        double y = s + sigma * randn(); /* noisy sample */

        /* LLR for soft-decision Viterbi */
        LLR[i] = (2.0 * y) / var;

        /* Hard decision */
        rx_bits[i] = (y >= 0.0 ? 0 : 1);
      }

      /* --- 4) Decode --- */
      nsc_decode_r05_soft(LLR, info_soft, code_hat);
      nsc_decode_r05_hard(rx_bits, info_hard, code_hat);

      /* --- 5) Count bit errors --- */
      for (int i = 0; i < K; i++) {
        if (info_soft[i] != data[i])
          error_soft++;
        if (info_hard[i] != data[i])
          error_hard++;
      }

      total_bits += K;
    }

    /* ------------------------------------------------------------------
     *  BER evaluation
     * ------------------------------------------------------------------ */
    double BER_soft = (double)error_soft / total_bits;
    double BER_hard = (double)error_hard / total_bits;
    double BER_bpsk = bpsk_ber(EbN0);

    printf("%.1f, %.10f, %.10f, %.10f\n", EbN0_dB, BER_soft, BER_hard,
           BER_bpsk);

    fprintf(fp, "%.1f,%.10f,%.10f,%.10f\n", EbN0_dB, BER_soft, BER_hard,
            BER_bpsk);
  }

  /* Cleanup */
  fclose(fp);

  free(data);
  free(code);
  free(LLR);
  free(rx_bits);
  free(info_soft);
  free(info_hard);
  free(code_hat);

  return 0;
}
