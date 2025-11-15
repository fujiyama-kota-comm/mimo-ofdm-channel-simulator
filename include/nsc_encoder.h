#ifndef NSC_ENCODER_H
#define NSC_ENCODER_H

#ifdef __cplusplus
extern "C" {
#endif

/* =============================================================================
 *  NSC Encoder (Non-Systematic Convolutional Code) — Rate 1/2
 * =============================================================================
 *
 *  Header for the table-driven NSC encoder operating on:
 *
 *      • Constraint length K = 3
 *      • 4 states (2-bit shift-register)
 *      • Non-systematic output
 *      • Rate 1/2 (1 input bit → 2 coded bits)
 *
 *  The encoder implementation (nsc_encoder.c) relies exclusively on the
 *  trellis tables defined in trellis.h:
 *
 *      - nsc_output_bits[state][input_bit][2]
 *      - nsc_next_state[state][input_bit]
 *
 *  No conditional branching is used in the main encoding loop.
 *
 *  ---------------------------------------------------------------------------
 *  Termination
 *  ---------------------------------------------------------------------------
 *  Two tail bits {0,0} are appended internally to force the trellis back to
 *  STATE_A (00). This simplifies Viterbi traceback and ensures a well-defined
 *  final state.
 *
 *  ---------------------------------------------------------------------------
 *  Output Length
 *  ---------------------------------------------------------------------------
 *      Input length  = K  (nsc_info_len)
 *      Tail length   = 2  (nsc_tail_len)
 *
 *      Output length = 2 * (K + 2)
 *
 *  The caller MUST set the following global parameters before encoding:
 *
 *      nsc_info_len   : number of input bits K
 *      nsc_code_len   : total number of output coded bits
 *      nsc_tail_len   : number of termination bits (default = 2)
 *
 * =============================================================================
 */

/* Global parameters (must be set by caller before using the encoder) --------
 */
extern int nsc_info_len; /* Number of information bits K */
extern int nsc_code_len; /* Output length = 2 * (K + nsc_tail_len) */
extern int nsc_tail_len; /* Tail bits (default: 2) */

/* =============================================================================
 *  Function: nsc_encode_r05
 * =============================================================================
 *  Encode input bits using a rate-1/2 non-systematic convolutional code.
 *
 *  Parameters:
 *      data : [input]  array of K information bits
 *      code : [output] array of length 2 * (K + nsc_tail_len)
 *
 *  Notes:
 *      - Two tail bits {0,0} are appended internally.
 *      - Encoding is fully table-driven (branchless).
 *
 * =============================================================================
 */
void nsc_encode_r05(const int *data, int *code);

#ifdef __cplusplus
}
#endif

#endif /* NSC_ENCODER_H */
