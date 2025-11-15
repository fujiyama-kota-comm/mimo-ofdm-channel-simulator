#include "trellis.h"

/* =============================================================================
 *  Trellis Definition for Non-Systematic Convolutional Code (NSC) — Rate 1/2
 * =============================================================================
 *
 *  This file defines the trellis structure for a 2-bit shift-register,
 *  rate-1/2, non-systematic convolutional encoder. The trellis is represented
 *  as lookup tables (no branching) for high-speed Viterbi decoding.
 *
 *  ---------------------------------------------------------------------------
 *  State Representation (2-bit shift register)
 *  ---------------------------------------------------------------------------
 *      STATE_A : 00
 *      STATE_B : 01
 *      STATE_C : 10
 *      STATE_D : 11
 *
 *  For input bit b ∈ {0,1}, the encoder outputs a 2-bit vector (v,w) and
 *  transitions to the next state, determined purely by the shift-register.
 *
 *  These tables enable:
 *      - Branchless encoder implementation
 *      - Fast metric computation in Viterbi decoding
 *      - Clear, deterministic trellis representation
 *
 *  ---------------------------------------------------------------------------
 *  NOTE:
 *      This trellis is NOT derived from arbitrary generator polynomials.
 *      Instead, it is a fixed “educational / demo” NSC example with:
 *
 *          - Constraint length K = 3 (2-bit state)
 *          - Non-systematic output (only (v,w), no systematic bit)
 *          - Two output bits per input bit → rate 1/2
 *
 *      This differs from conventional systematic encoders using polynomials
 *      such as 133, 171 (octal). Here, a simpler non-systematic form is used.
 *
 *  ---------------------------------------------------------------------------
 *  Provided lookup tables:
 *      nsc_output_bits[state][input][2] → output bits (v,w)
 *      nsc_next_state[state][input]     → next trellis state
 *
 *  These tables are consumed by both:
 *      - nsc_encoder.c
 *      - nsc_decoder.c (hard-decision & soft-decision Viterbi)
 *
 * =============================================================================
 */

/* =============================================================================
 *  Output Bit Table: nsc_output_bits[state][input_bit][2]
 * =============================================================================
 *  Example:
 *      nsc_output_bits[STATE_A][0] = {0,0}
 *      nsc_output_bits[STATE_A][1] = {1,1}
 *
 *  Meaning:
 *      At state 00, input 1 produces output 11 and moves to STATE_C (10).
 *
 *  Structure:
 *      state = 0..3 (A..D)
 *      input = 0 or 1
 *      output = {v, w} (two coded bits)
 *
 * =============================================================================
 */
const int nsc_output_bits[4][2][2] = {

    /* STATE_A : 00 ----------------------------------------------------------
     */
    {
        {0, 0}, /* input 0 → output 00 → next = STATE_A */
        {1, 1}  /* input 1 → output 11 → next = STATE_C */
    },

    /* STATE_B : 01 ----------------------------------------------------------
     */
    {
        {1, 1}, /* input 0 → output 11 → next = STATE_A */
        {0, 0}  /* input 1 → output 00 → next = STATE_C */
    },

    /* STATE_C : 10 ----------------------------------------------------------
     */
    {
        {0, 1}, /* input 0 → output 01 → next = STATE_B */
        {1, 0}  /* input 1 → output 10 → next = STATE_D */
    },

    /* STATE_D : 11 ----------------------------------------------------------
     */
    {
        {1, 0}, /* input 0 → output 10 → next = STATE_B */
        {0, 1}  /* input 1 → output 01 → next = STATE_D */
    }};

/* =============================================================================
 *  Next-State Table: nsc_next_state[state][input_bit]
 * =============================================================================
 *
 *  Shift register update rule:
 *
 *      new_state = ((state << 1) + input_bit) & 0b11
 *
 *  But for faster decoding, the next states are explicitly stored.
 *
 *  This table is essential for:
 *      - Encoder state tracking
 *      - Viterbi forward path computation
 *
 * =============================================================================
 */
const int nsc_next_state[4][2] = {
    /* STATE_A (00) */ {STATE_A, STATE_C},
    /* STATE_B (01) */ {STATE_A, STATE_C},
    /* STATE_C (10) */ {STATE_B, STATE_D},
    /* STATE_D (11) */ {STATE_B, STATE_D}};
