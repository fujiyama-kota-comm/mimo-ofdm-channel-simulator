#ifndef TRELLIS_H
#define TRELLIS_H

#ifdef __cplusplus
extern "C" {
#endif

/* =============================================================================
 *  Trellis Definition for Non-Systematic Convolutional Code (NSC) — Rate 1/2
 * =============================================================================
 *
 *  This header declares the trellis structure used by both the encoder and
 *  Viterbi decoder of a 2-bit shift-register, 4-state, non-systematic
 *  convolutional code.
 *
 *  The implementation is fully table-driven:
 *
 *      - NSCState                 : Enumeration of the 4 trellis states
 *      - nsc_output_bits[state][input_bit][2]
 *      - nsc_next_state[state][input_bit]
 *
 *  These tables are defined in trellis.c and are used by:
 *
 *      - nsc_encoder.c (branchless encoding)
 *      - nsc_decoder.c (hard/soft Viterbi decoding)
 *
 *  ---------------------------------------------------------------------------
 *  State Representation (2-bit shift register)
 *  ---------------------------------------------------------------------------
 *      STATE_A : 00
 *      STATE_B : 01
 *      STATE_C : 10
 *      STATE_D : 11
 *
 *  For each input bit b ∈ {0,1}, the trellis specifies:
 *
 *      (1) Output coded bits (v,w)
 *      (2) Next state of the shift register
 *
 *  The tables in this header provide a minimal and efficient interface for
 *  convolutional encoding and Viterbi metric computations.
 *
 * =============================================================================
 */

/* -----------------------------------------------------------------------------
 *  Enumeration of the 4 trellis states (2-bit shift-register contents)
 * ---------------------------------------------------------------------------
 */
typedef enum {
  STATE_A = 0, /* register = 00 */
  STATE_B = 1, /* register = 01 */
  STATE_C = 2, /* register = 10 */
  STATE_D = 3  /* register = 11 */
} NSCState;

/* =============================================================================
 *  Output Bit Lookup Table
 *      nsc_output_bits[state][input_bit][2]
 * =============================================================================
 *
 *  Output bits (v,w) produced by the convolutional encoder for each pair:
 *
 *      state ∈ {STATE_A, STATE_B, STATE_C, STATE_D}
 *      input_bit ∈ {0,1}
 *
 *  Access example:
 *      int v = nsc_output_bits[state][input][0];
 *      int w = nsc_output_bits[state][input][1];
 *
 *  This table is used by:
 *      - NSC encoder
 *      - Viterbi branch metric computation (hard / soft decision)
 *
 * =============================================================================
 */
extern const int nsc_output_bits[4][2][2];

/* =============================================================================
 *  Next-State Lookup Table
 *      nsc_next_state[state][input_bit]
 * =============================================================================
 *
 *  Determines the next state of the shift register:
 *
 *      new_state = nsc_next_state[current_state][input_bit]
 *
 *  This table eliminates all conditionals inside the encoder and Viterbi
 *  decoder, enabling efficient branchless implementations.
 *
 * =============================================================================
 */
extern const int nsc_next_state[4][2];

#ifdef __cplusplus
}
#endif

#endif /* TRELLIS_H */
