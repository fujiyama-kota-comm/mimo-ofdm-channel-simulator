/**
 * @file channel_tdl.c
 * @brief 3GPP TDL-A/B/C channel model with Rayleigh fading + Jakes-like
 *        time correlation (AR(1) using J0 autocorrelation).
 *
 * Profile:
 *   - delay_norm[i] : normalized delay τ_model,n (Table 7.7.2-1 "Normalized
 * delay")
 *   - power_lin[i]  : 10^(Power_dB/10)
 *
 * Usage:
 *   tdl_channel_t *ch = tdl_create(TDL_A, fs, fc, speed_kmh);
 *   ...
 *   tdl_update(ch);  // update fading taps h[l]
 *   ...
 *   tdl_get_channel(ch, delay_sec_array, cir, N);
 */

#include "channel_tdl.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846

/* ------------------------- Gaussian ------------------------- */
/* Zero-mean, unit-variance Gaussian using Box–Muller */
static float randn(void) {
  float u1 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 2.0f);
  float u2 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 2.0f);
  return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * PI * u2);
}

/* ------------------------- Doppler -------------------------- */
static double calc_fd(double fc, double speed_kmh) {
  double v = speed_kmh / 3.6; /* [m/s] */
  double lambda = 3e8 / fc;
  return v / lambda;
}

/* dB → linear power */
static double db2lin(double db) { return pow(10.0, db / 10.0); }

/* Bessel J0 (very simple approximation, sufficient for fd*Ts small) */
static double bessel_j0(double x) {
  double ax = fabs(x);
  if (ax < 8.0) {
    double y = x * x;
    return 1.0 - y * (0.25 - y * (0.046875 - y * 0.00390625));
  } else {
    double z = 8.0 / ax;
    double xx = ax - 0.7853981633974483; /* π/4 */
    return sqrt(0.636619772 / ax) * (cos(xx) - z * sin(xx));
  }
}

/* ------------------------------------------------------------
 * load_profile():
 *   Fill delay_norm[] and power_lin[] from 3GPP TDL-A/B/C tables.
 *
 *   delay_norm[i] : normalized delay τ_model,n (dimensionless)
 *   power_lin[i]  : linear power = 10^(Power_dB/10)
 * ------------------------------------------------------------ */
static void load_profile(tdl_channel_t *ch, tdl_type_t type) {
  int N = 0;

  if (type == TDL_A) {
    N = TDL_A_TAPS;
    ch->delay_norm = malloc(sizeof(double) * N);
    ch->power_lin = malloc(sizeof(double) * N);
    if (!ch->delay_norm || !ch->power_lin) {
      fprintf(stderr, "load_profile: malloc failed (TDL_A)\n");
      exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++) {
      ch->delay_norm[i] = TDL_A_DELAY[i]; /* τ_model,n (normalized) */
      ch->power_lin[i] = db2lin(TDL_A_POWER_DB[i]);
    }

  } else if (type == TDL_B) {
    N = TDL_B_TAPS;
    ch->delay_norm = malloc(sizeof(double) * N);
    ch->power_lin = malloc(sizeof(double) * N);
    if (!ch->delay_norm || !ch->power_lin) {
      fprintf(stderr, "load_profile: malloc failed (TDL_B)\n");
      exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++) {
      ch->delay_norm[i] = TDL_B_DELAY[i];
      ch->power_lin[i] = db2lin(TDL_B_POWER_DB[i]);
    }

  } else if (type == TDL_C) {
    N = TDL_C_TAPS;
    ch->delay_norm = malloc(sizeof(double) * N);
    ch->power_lin = malloc(sizeof(double) * N);
    if (!ch->delay_norm || !ch->power_lin) {
      fprintf(stderr, "load_profile: malloc failed (TDL_C)\n");
      exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++) {
      ch->delay_norm[i] = TDL_C_DELAY[i];
      ch->power_lin[i] = db2lin(TDL_C_POWER_DB[i]);
    }
  } else {
    fprintf(stderr, "load_profile: unsupported TDL type=%d\n", type);
    exit(EXIT_FAILURE);
  }

  ch->num_taps = N;
}

/* ------------------------------------------------------------
 * Create TDL channel object.
 *
 *   - delay_norm[] : normalized delays τ_model,n
 *   - power_lin[]  : linear powers
 *   - h[]          : fading taps (Rayleigh) with Jakes-like correlation
 * ------------------------------------------------------------ */
tdl_channel_t *tdl_create(tdl_type_t type, double fs, double fc,
                          double speed_kmh) {
  tdl_channel_t *ch = malloc(sizeof(tdl_channel_t));
  if (!ch) {
    fprintf(stderr, "tdl_create: malloc failed\n");
    return NULL;
  }

  memset(ch, 0, sizeof(tdl_channel_t));

  ch->type = type;
  ch->fs = fs;
  ch->fc = fc;
  ch->speed_kmh = speed_kmh;

  load_profile(ch, type);

  /* fading memory */
  ch->h = malloc(sizeof(complex float) * ch->num_taps);
  if (!ch->h) {
    fprintf(stderr, "tdl_create: malloc failed for h\n");
    free(ch->delay_norm);
    free(ch->power_lin);
    free(ch);
    return NULL;
  }

  /* Initialize taps as Rayleigh (CN(0, power_lin[i])) */
  for (int i = 0; i < ch->num_taps; i++) {
    float re = randn();
    float im = randn();
    float scale = sqrtf(ch->power_lin[i] * 0.5f);
    ch->h[i] = scale * (re + I * im);
  }
  return ch;
}

/* ------------------------------------------------------------
 * Free TDL channel object.
 * ------------------------------------------------------------ */
void tdl_free(tdl_channel_t *ch) {
  if (!ch)
    return;
  free(ch->delay_norm);
  free(ch->power_lin);
  free(ch->h);
  free(ch);
}

/* ------------------------------------------------------------
 * TDL Update (Rayleigh + Jakes-like time correlation)
 *
 *   h_i[n+1] = ρ h_i[n] + σ w_i[n],
 *   ρ = J0(2π fd Ts),  σ = sqrt(1 - ρ^2)
 *
 *   w_i[n] ~ CN(0, power_lin[i])
 * ------------------------------------------------------------ */
void tdl_update(tdl_channel_t *ch) {
  double dt = 1.0 / ch->fs;
  double fd = calc_fd(ch->fc, ch->speed_kmh);

  double rho = bessel_j0(2.0 * PI * fd * dt);
  double sigma = sqrt(fmax(0.0, 1.0 - rho * rho));

  for (int i = 0; i < ch->num_taps; i++) {
    float re = randn();
    float im = randn();
    float scale = sqrtf(ch->power_lin[i] * 0.5f);
    complex float w = scale * (re + I * im);

    ch->h[i] = rho * ch->h[i] + (float)sigma * w;
  }
}

/* ------------------------------------------------------------
 * CIR generator (optional helper)
 *
 * Build time-domain CIR out[n] from tap coefficients h[l] and
 * per-tap delays delay_sec[l] (in [s]).
 *
 * Example usage:
 *   double delay_sec[L];
 *   for (l) {
 *       double tau_ns = ch->delay_norm[l] * DS_desired_ns;
 *       delay_sec[l] = tau_ns * 1e-9;
 *   }
 *   tdl_get_channel(ch, delay_sec, cir, N);
 * ------------------------------------------------------------ */
void tdl_get_channel(const tdl_channel_t *ch,
                     const double *delay_sec, /* [s] per tap */
                     complex float *out, int N) {
  /* Initialize CIR */
  for (int n = 0; n < N; n++)
    out[n] = 0.0f + 0.0f * I;

  /* Place taps at rounded sample delays */
  for (int l = 0; l < ch->num_taps; l++) {
    int d = (int)lround(delay_sec[l] * ch->fs); /* seconds → samples */
    if (d >= 0 && d < N)
      out[d] += ch->h[l];
  }
}
