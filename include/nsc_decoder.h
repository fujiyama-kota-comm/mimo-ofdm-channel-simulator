#ifndef NSC_DECODER_H
#define NSC_DECODER_H

#ifdef __cplusplus
extern "C" {
#endif

/* =============================================================================
 *  NSC Decoder (Non-Systematic Convolutional Code) — Rate 1/2
 * =============================================================================
 *
 *  This header declares the Viterbi decoders for the non-systematic
 *  convolutional code defined in trellis.h.
 *
 *      • 4 trellis states (2-bit shift register)
 *      • Constraint length K = 3
 *      • Non-systematic rate-1/2 output
 *      • Terminated code: tail bits {0,0} force start/end state = STATE_A
 *
 *  Supported decoding modes:
 *
 *      (1) Soft-decision Viterbi : LLR input (double precision)
 *      (2) Hard-decision Viterbi : binary input (0/1 values)
 *
 *  Input/Output lengths (must be configured by caller):
 *
 *      nsc_info_len  = K
 *      nsc_code_len  = 2 * (K + nsc_tail_len)
 *
 *      Number of trellis steps = K + nsc_tail_len
 *
 *  During traceback, the last nsc_tail_len steps are discarded and only
 *  the first K decisions are written into info_hat[].
 *
 * =============================================================================
 */

/* =============================================================================
 *  Soft-Decision Viterbi Decoder
 * =============================================================================
 *
 *  Parameters:
 *      LLR      : input LLR array of length N (double)
 *      info_hat : output estimated information bits (length K)
 *      code_hat : optional re-encoded sequence (NULL → unused)
 *
 *  Notes:
 *      - Uses a negative correlation metric based on BPSK LLRs.
 *      - Performs standard Viterbi forward recursion + traceback.
 *
 * =============================================================================
 */
void nsc_decode_r05_soft(const double *LLR, int *info_hat, int *code_hat);

/* =============================================================================
 *  Hard-Decision Viterbi Decoder
 * =============================================================================
 *
 *  Parameters:
 *      rx_bits  : input hard-detected bits (0/1), length N
 *      info_hat : output estimated information bits (length K)
 *      code_hat : optional re-encoded sequence (NULL → unused)
 *
 *  Notes:
 *      - Uses Hamming distance (0,1,2) as the branch metric.
 *      - Otherwise identical to the soft-decision algorithm.
 *
 * =============================================================================
 */
void nsc_decode_r05_hard(const int *rx_bits, int *info_hat, int *code_hat);

#ifdef __cplusplus
}
#endif

#endif /* NSC_DECODER_H */
