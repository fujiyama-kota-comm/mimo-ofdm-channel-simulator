/**
 * @file mimo_ofdm_ber.c
 * @brief 2x2 MIMO-OFDM BER simulation over 3GPP TDL channel
 *        with Jakes-like time correlation and linear equalizers.
 *
 * Features:
 *  - Nfft = 256, SCS = 30 kHz → fs = Nfft * SCS
 *  - MIMO: 2 Tx × 2 Rx spatial multiplexing (single user)
 *  - Channel: per-link 3GPP TDL-A/B/C (independent Rayleigh fading)
 *  - Equalizers:
 *      - ZF :   X_hat = (H^H H)^-1 H^H Y
 *      - MMSE: X_hat = (H^H H + σ² I)^-1 H^H Y
 *  - Modulations: BPSK, QPSK, 16/64/256QAM (Gray, Es=1 per stream)
 *  - Perfect CSI: frequency-domain H[k] is built from true taps
 *  - Output: results/mimo_ofdm_ber_data.csv
 *      Columns: SNR_dB, BER_BPSK, BER_QPSK, BER_16QAM, ...
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

#include "channel_tdl.h"
#include "mimo.h"
#include "modulator.h"
#include "ofdm.h"

#define PI 3.14159265358979323846

/* ------------------------------------------------------------
 * MIMO & OFDM Parameters
 * ------------------------------------------------------------ */
const int Nfft = 256;
const float SCS_HZ = 30e3f;
const double fs = (double)Nfft * (double)SCS_HZ;

const int Ntx = 2; /* number of transmit antennas */
const int Nrx = 2; /* number of receive antennas */

const tdl_type_t tdl_type = TDL_A;
const double fc = 3.5e9;
const double speed_kmh = 30.0;

/* RMS delay spread scaling [ns] */
const double DS_desired_ns = 363.0;

/* Modulations to simulate */
modulation_t mod_list[] = {MOD_BPSK, MOD_QPSK, MOD_16QAM, MOD_64QAM,
                           MOD_256QAM};
const int NUM_MOD = sizeof(mod_list) / sizeof(mod_list[0]);

/* SNR sweep (Es/N0 in dB) */
const float SNR_MIN_dB = 0.0;
const float SNR_MAX_dB = 50.0;
const float SNR_STEP_dB = 5.0;

/* Monte Carlo iterations per SNR & modulation */
const int Niter = 100000;

/* Equalizer selection */
typedef enum { EQ_ZF = 0, EQ_MMSE = 1 } eq_method_t;
const eq_method_t EQ_METHOD = EQ_MMSE; /* change to EQ_ZF if needed */

/* ------------------------------------------------------------
 * Utility: Gaussian & AWGN
 * ------------------------------------------------------------ */
static float randn(void) {
  float u1 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 2.0f);
  float u2 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 2.0f);
  return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * PI * u2);
}

/* Add complex AWGN assuming Es = 1 per stream.
 * SNR_dB is Es/N0 in dB.
 */
static void add_awgn_freq(complex float *X, int N, float snr_db) {
  float snr = powf(10.0f, snr_db / 10.0f);
  float sigma = sqrtf(1.0f / (2.0f * snr)); /* per real dimension */
  for (int k = 0; k < N; k++) {
    X[k] += sigma * randn() + I * sigma * randn();
  }
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

/* String name for modulation (for logging / CSV header) */
static const char *modulation_name(modulation_t mod) {
  switch (mod) {
  case MOD_BPSK:
    return "BPSK";
  case MOD_QPSK:
    return "QPSK";
  case MOD_16QAM:
    return "16QAM";
  case MOD_64QAM:
    return "64QAM";
  case MOD_256QAM:
    return "256QAM";
  default:
    return "UNKNOWN";
  }
}

/* ------------------------------------------------------------
 * SISO-like TDL convolution: y[n] = Σ_l h[l] x[n - delay_samp[l]]
 * (single Tx-Rx link)
 * ------------------------------------------------------------ */
static void apply_tdl_fading_delay(const complex float *x, int Nsym,
                                   const complex float *h_tap,
                                   const int *delay_samp, int L,
                                   complex float *y) {
  for (int n = 0; n < Nsym; n++) {
    complex float sum = 0.0f + 0.0f * I;
    for (int l = 0; l < L; l++) {
      int idx = n - delay_samp[l];
      if (idx >= 0 && idx < Nsym) {
        sum += h_tap[l] * x[idx];
      }
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

  /* Open CSV */
  FILE *fp = fopen("results/mimo_ofdm_ber_data.csv", "w");
  if (!fp) {
    fprintf(stderr, "Cannot open results/mimo_ofdm_ber_data.csv\n");
    return EXIT_FAILURE;
  }

  /* CSV header */
  fprintf(fp, "SNR_dB");
  for (int mi = 0; mi < NUM_MOD; mi++) {
    fprintf(fp, ",BER_%s", modulation_name(mod_list[mi]));
  }
  fprintf(fp, "\n");

  /* Print basic info */
  printf("====================================================\n");
  printf(" MIMO-OFDM BER Simulation (2x2, TDL + Jakes + Perfect CSI)\n");
  printf("  Nfft = %d, SCS = %.1f kHz, fs = %.2f MHz\n", Nfft, SCS_HZ / 1e3,
         fs / 1e6);
  printf("  Channel = TDL-%c, Desired RMS DS = %.1f ns\n",
         (tdl_type == TDL_A)   ? 'A'
         : (tdl_type == TDL_B) ? 'B'
                               : 'C',
         DS_desired_ns);
  printf("  MIMO: %d Tx × %d Rx\n", Ntx, Nrx);
  printf("  Equalizer: %s\n", (EQ_METHOD == EQ_ZF) ? "ZF" : "MMSE");
  printf("====================================================\n");

  /* ----------------------------------------------------------
   * Create per-link TDL channels: H[rx][tx]
   * -------------------------------------------------------- */
  tdl_channel_t *tdl[Nrx][Ntx];

  for (int r = 0; r < Nrx; r++) {
    for (int t = 0; t < Ntx; t++) {
      tdl[r][t] = tdl_create(tdl_type, fs, fc, speed_kmh);
      if (!tdl[r][t]) {
        fprintf(stderr, "tdl_create failed (r=%d, t=%d)\n", r, t);
        return EXIT_FAILURE;
      }
    }
  }

  int L = tdl[0][0]->num_taps;

  /* Delay scaling (normalized τ_model → seconds → samples) */
  int *delay_samp = malloc(sizeof(int) * L);
  if (!delay_samp) {
    fprintf(stderr, "malloc failed for delay_samp\n");
    return EXIT_FAILURE;
  }

  int max_delay_samp = 0;
  for (int l = 0; l < L; l++) {
    double tau_model = tdl[0][0]->delay_norm[l];
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

  /* ----------------------------------------------------------
   * Allocate buffers
   * -------------------------------------------------------- */
  /* bits: [tx][bit_index] (max 8 bits/sym) */
  int **bits_in = malloc(Ntx * sizeof(int *));
  int **bits_out = malloc(Ntx * sizeof(int *));
  if (!bits_in || !bits_out) {
    fprintf(stderr, "malloc failed for bits arrays\n");
    return EXIT_FAILURE;
  }
  for (int t = 0; t < Ntx; t++) {
    bits_in[t] = malloc(sizeof(int) * Nfft * 8);
    bits_out[t] = malloc(sizeof(int) * Nfft * 8);
  }

  /* Frequency-domain symbols: X_tx[tx][k], X_rx[rx][k], X_eq_tx[tx][k] */
  complex float **X_tx = malloc(Ntx * sizeof(complex float *));
  complex float **X_eq_tx = malloc(Ntx * sizeof(complex float *));
  complex float **X_rx = malloc(Nrx * sizeof(complex float *));
  if (!X_tx || !X_eq_tx || !X_rx) {
    fprintf(stderr, "malloc failed for X_* arrays\n");
    return EXIT_FAILURE;
  }
  for (int t = 0; t < Ntx; t++) {
    X_tx[t] = malloc(sizeof(complex float) * Nfft);
    X_eq_tx[t] = malloc(sizeof(complex float) * Nfft);
  }
  for (int r = 0; r < Nrx; r++) {
    X_rx[r] = malloc(sizeof(complex float) * Nfft);
  }

  /* Time-domain symbols: x_time_tx[tx][n], y_time_rx[rx][n] */
  complex float **x_time_tx = malloc(Ntx * sizeof(complex float *));
  complex float **y_time_rx = malloc(Nrx * sizeof(complex float *));
  if (!x_time_tx || !y_time_rx) {
    fprintf(stderr, "malloc failed for time-domain arrays\n");
    return EXIT_FAILURE;
  }
  for (int t = 0; t < Ntx; t++) {
    x_time_tx[t] = malloc(sizeof(complex float) * Nsym_time);
  }
  for (int r = 0; r < Nrx; r++) {
    y_time_rx[r] = malloc(sizeof(complex float) * Nsym_time);
  }

  /* Temp buffer for single-link convolution */
  complex float *tmp_conv = malloc(sizeof(complex float) * Nsym_time);
  if (!tmp_conv) {
    fprintf(stderr, "malloc failed for tmp_conv\n");
    return EXIT_FAILURE;
  }

  /* For building per-link frequency responses H_rt[rx,tx,k] */
  complex float *h_time = malloc(sizeof(complex float) * Nfft);
  complex float *H_rt = malloc(sizeof(complex float) * Nrx * Ntx * Nfft);
  if (!h_time || !H_rt) {
    fprintf(stderr, "malloc failed for H_rt buffers\n");
    return EXIT_FAILURE;
  }
#define HRT_IDX(r, t, k) (((r) * Ntx + (t)) * Nfft + (k))

  /* MIMO equalizer interface: H_f[k], Y_f[k], X_hat[k] */
  mimo_matrix_t **H_f = malloc(Nfft * sizeof(mimo_matrix_t *));
  mimo_matrix_t **Y_f = malloc(Nfft * sizeof(mimo_matrix_t *));
  mimo_matrix_t **X_hat = malloc(Nfft * sizeof(mimo_matrix_t *));
  if (!H_f || !Y_f || !X_hat) {
    fprintf(stderr, "malloc failed for MIMO matrix arrays\n");
    return EXIT_FAILURE;
  }

  for (int k = 0; k < Nfft; k++) {
    H_f[k] = mimo_alloc_matrix(Nrx, Ntx);
    Y_f[k] = mimo_alloc_matrix(Nrx, 1);
    X_hat[k] = mimo_alloc_matrix(Ntx, 1);
  }

  /* ============================================================
   * LOOP over SNR
   * ============================================================ */
  for (float snr_db = SNR_MIN_dB; snr_db <= SNR_MAX_dB + 1e-6f;
       snr_db += SNR_STEP_dB) {

    double ber_mod[NUM_MOD];
    memset(ber_mod, 0, sizeof(ber_mod));

    printf("\n=== SNR = %.1f dB ===\n", snr_db);

    float snr_lin = powf(10.0f, snr_db / 10.0f);
    double noise_var = 1.0 / (double)snr_lin; /* σ² (Es=1 per stream) */

    /* ----------------------------------------------------------
     * LOOP over modulation
     * -------------------------------------------------------- */
    for (int mi = 0; mi < NUM_MOD; mi++) {

      modulation_t mod = mod_list[mi];
      int bps = bits_per_symbol(mod);
      int n_bits = bps * Nfft;

      long long err = 0;
      long long total = 0;

      printf("  %s ... ", modulation_name(mod));
      fflush(stdout);

      for (int it = 0; it < Niter; it++) {

        /* 1) Generate bits & map to frequency-domain symbols */
        for (int t = 0; t < Ntx; t++) {
          for (int i = 0; i < n_bits; i++) {
            bits_in[t][i] = rand() & 1;
          }
          mod_map(bits_in[t], X_tx[t], Nfft, mod);
        }

        /* 2) OFDM modulation for each Tx antenna */
        for (int t = 0; t < Ntx; t++) {
          ofdm_modulate_symbol(X_tx[t], x_time_tx[t], &ofdm_cfg);
        }

        /* 3) MIMO TDL channel (time-domain convolution) */
        for (int r = 0; r < Nrx; r++) {
          for (int n = 0; n < Nsym_time; n++) {
            y_time_rx[r][n] = 0.0f + 0.0f * I;
          }
        }

        /* Update all TDL taps and propagate */
        for (int r = 0; r < Nrx; r++) {
          for (int t = 0; t < Ntx; t++) {
            tdl_update(tdl[r][t]);
            apply_tdl_fading_delay(x_time_tx[t], Nsym_time, tdl[r][t]->h,
                                   delay_samp, L, tmp_conv);
            for (int n = 0; n < Nsym_time; n++) {
              y_time_rx[r][n] += tmp_conv[n];
            }
          }
        }

        /* 4) OFDM demodulation at each Rx antenna */
        for (int r = 0; r < Nrx; r++) {
          ofdm_demodulate_symbol(y_time_rx[r], X_rx[r], &ofdm_cfg);
        }

        /* 5) Add AWGN */
        for (int r = 0; r < Nrx; r++) {
          add_awgn_freq(X_rx[r], Nfft, snr_db);
        }

        /* 6) Build per-link frequency responses H_rt[rx,tx,k] */
        for (int r = 0; r < Nrx; r++) {
          for (int t = 0; t < Ntx; t++) {
            for (int n = 0; n < Nfft; n++) {
              h_time[n] = 0.0f + 0.0f * I;
            }
            for (int l = 0; l < L; l++) {
              int d = delay_samp[l];
              if (d < Nfft) {
                h_time[d] += tdl[r][t]->h[l];
              }
            }
            ofdm_fft_symbol(h_time, &H_rt[HRT_IDX(r, t, 0)], &ofdm_cfg);
          }
        }

        /* 7) Build per-subcarrier MIMO matrices H_f[k], Y_f[k] */
        for (int k = 0; k < Nfft; k++) {
          /* Y: Nrx x 1 */
          for (int r = 0; r < Nrx; r++) {
            Y_f[k]->data[r] = X_rx[r][k];
          }
          /* H: Nrx x Ntx */
          for (int r = 0; r < Nrx; r++) {
            for (int t = 0; t < Ntx; t++) {
              H_f[k]->data[r * Ntx + t] = H_rt[HRT_IDX(r, t, k)];
            }
          }
        }

        /* 8) MIMO equalization per subcarrier */
        if (EQ_METHOD == EQ_ZF) {
          mimo_equalize_zf(H_f, Y_f, X_hat, Nfft, Nrx, Ntx);
        } else {
          mimo_equalize_mmse(H_f, Y_f, X_hat, Nfft, Nrx, Ntx, noise_var, 1.0);
        }

        /* 9) Collect equalized symbols per Tx antenna */
        for (int t = 0; t < Ntx; t++) {
          for (int k = 0; k < Nfft; k++) {
            X_eq_tx[t][k] = X_hat[k]->data[t]; /* (t,0) entry */
          }
        }

        /* 10) Demap and count bit errors over all Tx streams */
        for (int t = 0; t < Ntx; t++) {
          mod_demod(X_eq_tx[t], bits_out[t], Nfft, mod);
          for (int i = 0; i < n_bits; i++) {
            if (bits_in[t][i] != bits_out[t][i]) {
              err++;
            }
          }
          total += n_bits;
        }
      } /* end Monte Carlo loop */

      ber_mod[mi] = (double)err / (double)total;
      printf("BER = %.3e\n", ber_mod[mi]);
    }

    /* ----------------------------------------------------------
     * Write one row to CSV
     * -------------------------------------------------------- */
    fprintf(fp, "%.1f", snr_db);
    for (int mi = 0; mi < NUM_MOD; mi++) {
      fprintf(fp, ",%.8e", ber_mod[mi]);
    }
    fprintf(fp, "\n");
  } /* end SNR loop */

  fclose(fp);

  /* Cleanup */
  for (int t = 0; t < Ntx; t++) {
    free(bits_in[t]);
    free(bits_out[t]);
  }
  free(bits_in);
  free(bits_out);

  for (int t = 0; t < Ntx; t++) {
    free(X_tx[t]);
    free(X_eq_tx[t]);
  }
  for (int r = 0; r < Nrx; r++) {
    free(X_rx[r]);
  }
  free(X_tx);
  free(X_eq_tx);
  free(X_rx);

  for (int t = 0; t < Ntx; t++) {
    free(x_time_tx[t]);
  }
  for (int r = 0; r < Nrx; r++) {
    free(y_time_rx[r]);
  }
  free(x_time_tx);
  free(y_time_rx);

  free(tmp_conv);
  free(h_time);
  free(H_rt);
  free(delay_samp);

  for (int k = 0; k < Nfft; k++) {
    mimo_free_matrix(H_f[k]);
    mimo_free_matrix(Y_f[k]);
    mimo_free_matrix(X_hat[k]);
  }
  free(H_f);
  free(Y_f);
  free(X_hat);

  for (int r = 0; r < Nrx; r++) {
    for (int t = 0; t < Ntx; t++) {
      tdl_free(tdl[r][t]);
    }
  }

  return EXIT_SUCCESS;
}
