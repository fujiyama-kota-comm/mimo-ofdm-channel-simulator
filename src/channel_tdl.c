#include "channel_tdl.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846

/* ------------------------- Gaussian ------------------------- */
static float randn(void) {
  float u1 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 2.0f);
  float u2 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 2.0f);
  return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * PI * u2);
}

/* ------------------------- Doppler -------------------------- */
static double calc_fd(double fc, double speed_kmh) {
  double v = speed_kmh / 3.6; /* m/s */
  double lambda = 3e8 / fc;
  return v / lambda;
}

/* dB → linear */
static double db2lin(double db) { return pow(10.0, db / 10.0); }

/* Bessel J0 approximation */
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
 * load_profile(): 正規化遅延 τ_model,n とパワーを読み込む
 *   - delay_norm[i] : 38.901 Table 7.7.2-1 の "Normalized delay" そのまま
 *   - power_lin[i]  : 10^(dB/10)
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
      ch->delay_norm[i] = TDL_A_DELAY[i]; /* 無次元の τ_model,n */
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
  }

  ch->num_taps = N;
}

/* ------------------------------------------------------------
 * create TDL channel
 *   - delay_norm[] : 正規化遅延 τ_model,n
 *   - power_lin[]  : 線形パワー
 *   - h[]          : Rayleigh + Jakes 用のフェージングメモリ
 * ------------------------------------------------------------ */
tdl_channel_t *tdl_create(tdl_type_t type, double fs, double fc,
                          double speed_kmh) {
  tdl_channel_t *ch = malloc(sizeof(tdl_channel_t));
  if (!ch) {
    fprintf(stderr, "tdl_create: malloc failed for ch\n");
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

  for (int i = 0; i < ch->num_taps; i++) {
    float re = randn();
    float im = randn();
    float scale = sqrtf(ch->power_lin[i] * 0.5f);
    ch->h[i] = scale * (re + I * im);
  }
  return ch;
}

/* ------------------------------------------------------------
 * free
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
 * TDL Update (Rayleigh + Jakes)
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

    ch->h[i] = rho * ch->h[i] + sigma * w;
  }
}

/* ------------------------------------------------------------
 * CIR generator（オプション使用）
 *
 * ここでは「すでに秒 [s] に変換済みの delay_sec[]」を受け取り、
 * それをサンプルインデックスにして CIR を out[] に構成する。
 *
 * 例：
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
  for (int n = 0; n < N; n++)
    out[n] = 0.0f + 0.0f * I;

  for (int l = 0; l < ch->num_taps; l++) {
    int d = (int)lround(delay_sec[l] * ch->fs); /* 秒→サンプル */
    if (d >= 0 && d < N)
      out[d] += ch->h[l];
  }
}
