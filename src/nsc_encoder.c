#include "nsc_encoder.h"
#include "trellis.h"
#include <stdio.h>
#include <stdlib.h>

/* =============================================================================
 *  Non-Systematic Convolutional (NSC) Encoder — Rate 1/2
 * =============================================================================
 *
 *  This file implements a branchless, table-driven convolutional encoder
 *  operating on a 2-bit shift-register (4 states). The encoder depends
 *  solely on the trellis tables defined in trellis.h:
 *
 *      - nsc_output_bits[state][input_bit][2] : coded output bits (v,w)
 *      - nsc_next_state [state][input_bit]    : next shift-register state
 *
 *  The design is intentionally simple and efficient:
 *
 *      • Constraint length K = 3
 *      • Number of states   = 4 (00, 01, 10, 11)
 *      • Non-systematic output (no raw input bit in the output)
 *      • Rate 1/2 encoding (1 input bit → 2 output bits)
 *
 *  ---------------------------------------------------------------------------
 *  Tail-Biting / Termination
 *  ---------------------------------------------------------------------------
 *  Two tail bits {0,0} are appended to force the encoder back to the all-zero
 *  state (STATE_A). This is a common termination method for short-frame NSC
 *  codes and ensures compatibility with the Viterbi decoder's traceback.
 *
 *  ---------------------------------------------------------------------------
 *  Output Length
 *  ---------------------------------------------------------------------------
 *      input length  = K  (nsc_info_len)
 *      tail length   = 2  (nsc_tail_len)
 *
 *      total input bits  = K + 2
 *      total output bits = 2 * (K + 2)   // rate 1/2
 *
 *  IMPORTANT:
 *      The caller is responsible for setting:
 *
 *          nsc_info_len      (input length)
 *          nsc_code_len      (output length)
 *
 *      before calling the encoder.
 *
 * =============================================================================
 */

/* Global parameters (declared in header) ------------------------------------
 */
int nsc_info_len;     /* Number of information bits */
int nsc_code_len;     /* Output length = 2 * (nsc_info_len + nsc_tail_len) */
int nsc_tail_len = 2; /* Two termination bits forcing state → STATE_A */

/* =============================================================================
 *  nsc_encode_r05()
 * =============================================================================
 *
 *  Encode a binary sequence using a rate-1/2 NSC encoder.
 *
 *  Arguments:
 *      data : [input]  array of length K
 *      code : [output] array of length 2 * (K + 2)
 *
 *  Process:
 *      1) Copy input bits into a temporary buffer
 *      2) Append 2 tail bits = {0,0}
 *      3) Starting from STATE_A, follow the trellis transitions:
 *
 *              code[2*i]   = nsc_output_bits[state][b][0];
 *              code[2*i+1] = nsc_output_bits[state][b][1];
 *              state       = nsc_next_state[state][b];
 *
 *      This is a pure table-driven implementation: no branching, no bitwise
 *      logic required inside the main loop.
 *
 *  NOTE:
 *      Input size checks are intentionally omitted for performance.
 *
 * =============================================================================
 */
void nsc_encode_r05(const int *data, int *code) {
  int total_bits = nsc_info_len + nsc_tail_len; /* K + 2 */

  /* Allocate buffer for input + tail bits -------------------------------- */
  int *buffer = malloc(sizeof(int) * total_bits);
  if (!buffer) {
    fprintf(stderr, "[NSC Encoder] Memory allocation failed\n");
    return;
  }

  /* Copy information bits ------------------------------------------------- */
  for (int i = 0; i < nsc_info_len; i++) {
    buffer[i] = data[i];
  }

  /* Append 2 tail bits {0,0} to force encoder → STATE_A ------------------- */
  buffer[nsc_info_len] = 0;
  buffer[nsc_info_len + 1] = 0;

  /* Initial encoder state is always STATE_A (00) -------------------------- */
  int state = STATE_A;

  /* Main trellis loop: branchless encoding -------------------------------- */
  for (int i = 0; i < total_bits; i++) {

    int b = buffer[i]; /* input bit (0 or 1) */

    /* Write 2 coded output bits */
    code[2 * i] = nsc_output_bits[state][b][0];
    code[2 * i + 1] = nsc_output_bits[state][b][1];

    /* Move to next encoder state */
    state = nsc_next_state[state][b];
  }

  free(buffer);
}
