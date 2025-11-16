#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

#include "channel_tdl.h"
#include "modulator.h"
#include "ofdm.h"

#define PI 3.14159265358979323846

/* ------------------------------------------------------------
 * Simulation Parameters
 * ------------------------------------------------------------ */
const int Nfft = 256;
const float SCS_HZ = 30e3f;
const double fs = (double)Nfft * (double)SCS_HZ;

const tdl_type_t tdl_type = TDL_A;
const double fc = 3.5e9;
const double speed_kmh = 30.0;

/* RMS delay spread */
const double DS_desired_ns = 363.0;

/* Modulations */
modulation_t mod_list[] = {MOD_BPSK, MOD_QPSK, MOD_16QAM, MOD_64QAM,
                           MOD_256QAM};
const int NUM_MOD = sizeof(mod_list) / sizeof(mod_list[0]);

/* SNR sweep */
float SNRs[] = {0, 5, 10, 15, 20, 25, 30};
const int NUM_SNR = sizeof(SNRs) / sizeof(SNRs[0]);

/* Monte Carlo loops */
const int Niter = 3000;

/* Gaussian */
static float randn(void) {
  float u1 = ((float)rand() + 1) / ((float)RAND_MAX + 2);
  float u2 = ((float)rand() + 1) / ((float)RAND_MAX + 2);
  return sqrtf(-2 * logf(u1)) * cosf(2 * PI * u2);
}

/* AWGN */
static void add_awgn_freq(complex float *X, int N, float snr_db) {
  float snr = powf(10.0f, snr_db / 10.0f);
  float sigma = sqrtf(1.0f / (2.0f * snr));
  for (int k = 0; k < N; k++)
    X[k] += sigma * randn() + I * sigma * randn();
}

/* bits per symbol */
static int bits_per_symbol(modulation_t mod) {
  switch (mod) {
  case MOD_BPSK:
    return 1;
  case MOD_QPSK:
    return 2;
  case MOD_16QAM:
    return 4;
  case MOD_64QAM:
    return 6;
  case MOD_256QAM:
    return 8;
  default:
    return 0;
  }
}

/* TDL convolution */
static void apply_tdl_fading_delay(const complex float *x, int Nsym,
                                   const complex float *h_tap,
                                   const int *delay_samp, int L,
                                   complex float *y) {
  for (int n = 0; n < Nsym; n++) {
    complex float sum = 0;
    for (int l = 0; l < L; l++) {
      int idx = n - delay_samp[l];
      if (idx >= 0 && idx < Nsym)
        sum += h_tap[l] * x[idx];
    }
    y[n] = sum;
  }
}

/* ============================================================
 * MAIN
 * ============================================================ */
int main(void) {
  srand((unsigned)time(NULL));

#ifdef _WIN32
  _mkdir("results");
#else
  mkdir("results", 0777);
#endif

  /* Open CSV (horizontal format) */
  FILE *fp = fopen("results/ofdm_ber_data.csv", "w");
  if (!fp) {
    fprintf(stderr, "Cannot open results/ofdm_ber_data.csv\n");
    return -1;
  }

  /* CSV header */
  fprintf(fp, "SNR_dB");
  for (int mi = 0; mi < NUM_MOD; mi++) {
    fprintf(fp, ",BER_%dQAM", (int)mod_list[mi]);
  }
  fprintf(fp, "\n");

  /* Print info */
  printf("====================================================\n");
  printf(" OFDM BER Simulation (TDL + Jakes + Perfect CSI)\n");
  printf("  Nfft=%d, SCS=%.1f kHz, fs=%.2f MHz\n", Nfft, SCS_HZ / 1e3,
         fs / 1e6);
  printf("  Desired DS = %.1f ns\n", DS_desired_ns);
  printf("====================================================\n");

  tdl_channel_t *tdl = tdl_create(tdl_type, fs, fc, speed_kmh);
  int L = tdl->num_taps;

  /* Delay scaling */
  int *delay_samp = malloc(sizeof(int) * L);
  int max_delay_samp = 0;

  for (int l = 0; l < L; l++) {
    double tau_model = tdl->delay_norm[l];
    double tau_scaled_ns = tau_model * DS_desired_ns;
    double delay_sec = tau_scaled_ns * 1e-9;
    int d = (int)lround(delay_sec * fs);
    if (d < 0)
      d = 0;
    delay_samp[l] = d;
    if (d > max_delay_samp)
      max_delay_samp = d;
  }

  int Ncp = max_delay_samp + 4;
  int Nsym_time = Nfft + Ncp;

  ofdm_cfg_t ofdm_cfg = {Nfft, Ncp};

  /* Buffers */
  int *bits_in = malloc(sizeof(int) * Nfft * 8);
  int *bits_out = malloc(sizeof(int) * Nfft * 8);
  complex float *X_tx = malloc(sizeof(complex float) * Nfft);
  complex float *X_rx = malloc(sizeof(complex float) * Nfft);
  complex float *X_eq = malloc(sizeof(complex float) * Nfft);
  complex float *x_time = malloc(sizeof(complex float) * Nsym_time);
  complex float *y_time = malloc(sizeof(complex float) * Nsym_time);
  complex float *h_tap = malloc(sizeof(complex float) * L);
  complex float *h_time = malloc(sizeof(complex float) * Nfft);
  complex float *H_f = malloc(sizeof(complex float) * Nfft);

  /* ============================================================
   * LOOP over SNR  (→ CSV の 1 行になる)
   * ============================================================ */
  for (int si = 0; si < NUM_SNR; si++) {

    float snr_db = SNRs[si];
    double ber_mod[NUM_MOD];
    memset(ber_mod, 0, sizeof(ber_mod)); /* ゼロ初期化 */

    printf("\n=== SNR = %.1f dB ===\n", snr_db);

    /* ========================================================
     * LOOP over modulation
     * ====================================================== */
    for (int mi = 0; mi < NUM_MOD; mi++) {

      modulation_t mod = mod_list[mi];
      int bps = bits_per_symbol(mod);
      int n_bits = bps * Nfft;
      long long err = 0, total = 0;

      printf("  Mod %dQAM... ", (int)mod);

      for (int it = 0; it < Niter; it++) {

        for (int i = 0; i < n_bits; i++)
          bits_in[i] = rand() & 1;

        mod_map(bits_in, X_tx, Nfft, mod);
        ofdm_modulate_symbol(X_tx, x_time, &ofdm_cfg);

        tdl_update(tdl);
        for (int l = 0; l < L; l++)
          h_tap[l] = tdl->h[l];

        apply_tdl_fading_delay(x_time, Nsym_time, h_tap, delay_samp, L, y_time);

        ofdm_demodulate_symbol(y_time, X_rx, &ofdm_cfg);
        add_awgn_freq(X_rx, Nfft, snr_db);

        for (int i = 0; i < Nfft; i++)
          h_time[i] = 0;

        for (int l = 0; l < L; l++) {
          int d = delay_samp[l];
          if (d < Nfft)
            h_time[d] += h_tap[l];
        }

        ofdm_fft_symbol(h_time, H_f, &ofdm_cfg);

        for (int k = 0; k < Nfft; k++) {
          float mag2 = crealf(H_f[k]) * crealf(H_f[k]) +
                       cimagf(H_f[k]) * cimagf(H_f[k]) + 1e-12f;
          X_eq[k] = X_rx[k] * conjf(H_f[k]) / mag2;
        }

        mod_demod(X_eq, bits_out, Nfft, mod);

        for (int i = 0; i < n_bits; i++)
          if (bits_in[i] != bits_out[i])
            err++;

        total += n_bits;
      }

      ber_mod[mi] = (double)err / total;
      printf("BER = %.3e\n", ber_mod[mi]);
    }

    /* ========================================================
     * CSV 1 行出力
     * ====================================================== */
    fprintf(fp, "%.1f", snr_db);
    for (int mi = 0; mi < NUM_MOD; mi++) {
      fprintf(fp, ",%.8e", ber_mod[mi]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  free(bits_in);
  free(bits_out);
  free(X_tx);
  free(X_rx);
  free(X_eq);
  free(x_time);
  free(y_time);
  free(h_tap);
  free(h_time);
  free(H_f);
  free(delay_samp);
  tdl_free(tdl);

  return 0;
}
