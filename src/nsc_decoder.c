#include "nsc_decoder.h"
#include "nsc_encoder.h"
#include "trellis.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* =============================================================================
 *  NSC Decoder (Viterbi, Rate 1/2, Terminated)
 * =============================================================================
 *
 *  This file implements both soft-decision and hard-decision Viterbi decoders
 *  for the non-systematic convolutional code (NSC) defined in trellis.h.
 *
 *      • 2-bit shift register → 4 states (STATE_A .. STATE_D)
 *      • Constraint length K = 3
 *      • Non-systematic rate-1/2 output (2 coded bits per input bit)
 *      • Terminated code: tail bits {0,0} force final state = STATE_A
 *
 *  The decoder follows the standard Viterbi procedure:
 *
 *      (1) Forward recursion (branch metric + path metric update)
 *      (2) Backward traceback (recover input bits)
 *      (3) Optional re-encoding for consistency checking
 *
 *  Path/state transitions use only the lookup tables:
 *
 *      nsc_next_state[s][b]
 *      nsc_output_bits[s][b][2]
 *
 *  No branching is used for encoder-side transitions; all logic is
 * table-driven.
 *
 * =============================================================================
 */

extern int nsc_info_len; /* Number of information bits K          */
extern int nsc_code_len; /* Number of received bits N = 2*(K+2)   */

/* =============================================================================
 *  Branch Metric (Soft Decision)
 * =============================================================================
 *
 *  LLR[] contains one LLR value per coded bit (BPSK):
 *
 *      LLR[k] ≈ log( p(y_k | +1) / p(y_k | -1) )
 *
 *  Mapping:
 *      bit 0 → BPSK +1
 *      bit 1 → BPSK -1
 *
 *  For the 2 output bits (v,w), the soft branch metric is defined as:
 *
 *      bm = -(s0 * LLR_v + s1 * LLR_w)
 *
 *  Lower metric → more likely path.
 *
 * =============================================================================
 */
static inline double branch_metric_soft_symbol(const double *LLR, int sym,
                                               int out0, int out1) {
  double v = LLR[2 * sym];
  double w = LLR[2 * sym + 1];

  double s0 = (out0 == 0) ? +1.0 : -1.0;
  double s1 = (out1 == 0) ? +1.0 : -1.0;

  return -(s0 * v + s1 * w);
}

/* =============================================================================
 *  Branch Metric (Hard Decision)
 * =============================================================================
 *
 *  rx_bits[] contains 0/1 hard-detected values.
 *  The branch metric is the Hamming distance between:
 *
 *      (received bits)  vs  (expected output bits)
 *
 * =============================================================================
 */
static inline int branch_metric_hard_symbol(const int *rx_bits, int sym,
                                            int out0, int out1) {
  int v = rx_bits[2 * sym];
  int w = rx_bits[2 * sym + 1];
  return (v != out0) + (w != out1);
}

/* =============================================================================
 *  Soft-Decision Viterbi Decoding
 * =============================================================================
 *
 *  Parameters:
 *      LLR      : input LLRs (length N)
 *      info_hat : output estimated information bits (length K)
 *      code_hat : optional; if non-NULL, encoder is invoked for re-encoding
 *
 *  Notes:
 *      - Trellis start and end states are assumed to be STATE_A due to tail
 * bits
 *      - Traceback ignores the last two steps (tail region)
 *
 * =============================================================================
 */
void nsc_decode_r05_soft(const double *LLR, int *info_hat, int *code_hat) {
  int K = nsc_info_len;
  int N = nsc_code_len;
  int steps = N / 2; /* K + tail_len */

  double metric_prev[4], metric_curr[4];

  int *prev_state = malloc(sizeof(int) * steps * 4);
  int *prev_bit = malloc(sizeof(int) * steps * 4);

  if (!prev_state || !prev_bit) {
    fprintf(stderr, "[NSC Decoder] Memory allocation failed (soft)\n");
    free(prev_state);
    free(prev_bit);
    return;
  }

  /* -------------------------------------------------------------------------
   * Initialization
   * ---------------------------------------------------------------------- */
  for (int s = 0; s < 4; s++)
    metric_prev[s] = 1e30;

  metric_prev[STATE_A] = 0.0; /* start in all-zero state */

  /* -------------------------------------------------------------------------
   * Forward recursion
   * ---------------------------------------------------------------------- */
  for (int i = 0; i < steps; i++) {

    for (int s = 0; s < 4; s++)
      metric_curr[s] = 1e30;

    for (int ps = 0; ps < 4; ps++) {
      if (metric_prev[ps] >= 1e29)
        continue; /* unreachable state */

      for (int b = 0; b < 2; b++) {

        int ns = nsc_next_state[ps][b];
        int o0 = nsc_output_bits[ps][b][0];
        int o1 = nsc_output_bits[ps][b][1];

        double bm = branch_metric_soft_symbol(LLR, i, o0, o1);
        double cand = metric_prev[ps] + bm;

        if (cand < metric_curr[ns]) {
          metric_curr[ns] = cand;
          prev_state[i * 4 + ns] = ps;
          prev_bit[i * 4 + ns] = b;
        }
      }
    }

    for (int s = 0; s < 4; s++)
      metric_prev[s] = metric_curr[s];
  }

  /* -------------------------------------------------------------------------
   * Final state selection (termination enforces STATE_A, but evaluate anyway)
   * ---------------------------------------------------------------------- */
  int state = STATE_A;
  double best = metric_prev[state];

  for (int s = 0; s < 4; s++) {
    if (metric_prev[s] < best) {
      best = metric_prev[s];
      state = s;
    }
  }

  /* -------------------------------------------------------------------------
   * Backward Traceback
   * ---------------------------------------------------------------------- */
  for (int i = steps - 1; i >= 0; i--) {
    int b = prev_bit[i * 4 + state];
    int ps = prev_state[i * 4 + state];

    if (i < K)
      info_hat[i] = b;

    state = ps;
  }

  /* Optional re-encoding */
  if (code_hat)
    nsc_encode_r05(info_hat, code_hat);

  free(prev_state);
  free(prev_bit);
}

/* =============================================================================
 *  Hard-Decision Viterbi Decoding
 * =============================================================================
 *
 *  Parameters:
 *      rx_bits  : hard-detected 0/1 bits (length N)
 *      info_hat : estimated information bits (length K)
 *      code_hat : optional; output of re-encoding (NULL → unused)
 *
 *  Same structure as the soft-decision version, but uses Hamming distance
 *  as the branch metric.
 *
 * =============================================================================
 */
void nsc_decode_r05_hard(const int *rx_bits, int *info_hat, int *code_hat) {
  int K = nsc_info_len;
  int N = nsc_code_len;
  int steps = N / 2;

  int metric_prev[4], metric_curr[4];

  int *prev_state = malloc(sizeof(int) * steps * 4);
  int *prev_bit = malloc(sizeof(int) * steps * 4);

  if (!prev_state || !prev_bit) {
    fprintf(stderr, "[NSC Decoder] Memory allocation failed (hard)\n");
    free(prev_state);
    free(prev_bit);
    return;
  }

  /* Initialization */
  for (int s = 0; s < 4; s++)
    metric_prev[s] = 1000000000;

  metric_prev[STATE_A] = 0;

  /* Forward recursion */
  for (int i = 0; i < steps; i++) {

    for (int s = 0; s < 4; s++)
      metric_curr[s] = 1000000000;

    for (int ps = 0; ps < 4; ps++) {
      if (metric_prev[ps] >= 100000000)
        continue;

      for (int b = 0; b < 2; b++) {

        int ns = nsc_next_state[ps][b];
        int o0 = nsc_output_bits[ps][b][0];
        int o1 = nsc_output_bits[ps][b][1];

        int bm = branch_metric_hard_symbol(rx_bits, i, o0, o1);
        int cand = metric_prev[ps] + bm;

        if (cand < metric_curr[ns]) {
          metric_curr[ns] = cand;
          prev_state[i * 4 + ns] = ps;
          prev_bit[i * 4 + ns] = b;
        }
      }
    }

    for (int s = 0; s < 4; s++)
      metric_prev[s] = metric_curr[s];
  }

  /* Final state selection */
  int state = STATE_A;
  int best = metric_prev[state];

  for (int s = 0; s < 4; s++) {
    if (metric_prev[s] < best) {
      best = metric_prev[s];
      state = s;
    }
  }

  /* Backward traceback */
  for (int i = steps - 1; i >= 0; i--) {
    int b = prev_bit[i * 4 + state];
    int ps = prev_state[i * 4 + state];

    if (i < K)
      info_hat[i] = b;

    state = ps;
  }

  if (code_hat)
    nsc_encode_r05(info_hat, code_hat);

  free(prev_state);
  free(prev_bit);
}
