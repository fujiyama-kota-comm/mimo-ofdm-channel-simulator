/**
 * @file channel_tdl.h
 * @brief 3GPP 38.901 TDL-A/B/C channel model
 *
 * This header defines:
 *   - Tap profiles (delay_norm[], power_dB[])
 *   - tdl_channel_t structure
 *   - Rayleigh fading with Jakes-like temporal correlation
 *
 * Note:
 *   delay_norm[i] is the *normalized delay* τ_model,n given in 3GPP 38.901.
 *   Actual time delay is computed by the caller as:
 *
 *      τ_scaled[i] = delay_norm[i] * DS_desired   (ns)
 *      delay_sec[i] = τ_scaled[i] * 1e-9
 *
 * The library itself does NOT assume any RMS delay; scaling must be done
 * externally to support arbitrary delay spreads.
 */

#ifndef CHANNEL_TDL_H
#define CHANNEL_TDL_H

#include <complex.h>

/* ------------------------------------------------------------
 * TDL Type Selector (3GPP 38.901)
 * ------------------------------------------------------------ */
typedef enum {
  TDL_A = 0, /**< Urban Macro (UMa) NLOS-like profile */
  TDL_B = 1, /**< UMi NLOS-like profile */
  TDL_C = 2  /**< RMa LOS-like profile */
} tdl_type_t;

/* ------------------------------------------------------------
 * 3GPP 38.901 TDL tap profiles
 * delay_norm[i] : normalized delay (dimensionless)
 * power_dB[i]   : power level in dB
 *
 * These are fixed profiles from Table 7.7.2-1.
 * ------------------------------------------------------------ */

/* -------------------- TDL-A -------------------- */
#define TDL_A_TAPS 23
static const double TDL_A_DELAY[TDL_A_TAPS] = {
    0.0,    0.3819, 0.4025, 0.5868, 0.4610, 0.5375, 0.6708, 0.5750,
    0.7618, 1.5375, 1.8978, 2.2242, 2.1718, 2.4942, 2.5119, 3.0582,
    4.0810, 4.4579, 4.5695, 4.7966, 5.0066, 5.3043, 9.6586};
static const double TDL_A_POWER_DB[TDL_A_TAPS] = {
    -13.4, 0.0,   -2.2,  -4.0,  -6.0,  -8.2,  -9.9,  -10.5,
    -7.5,  -15.9, -6.6,  -16.7, -12.4, -15.2, -10.8, -11.3,
    -12.7, -16.2, -18.3, -18.9, -16.6, -19.9, -29.7};

/* -------------------- TDL-B -------------------- */
#define TDL_B_TAPS 23
static const double TDL_B_DELAY[TDL_B_TAPS] = {
    0.0000, 0.1072, 0.2155, 0.2095, 0.2870, 0.2986, 0.3752, 0.5055,
    0.3681, 0.3697, 0.5700, 0.5283, 1.1021, 1.2756, 1.5474, 1.7842,
    2.0169, 2.8294, 3.0219, 3.6187, 4.1067, 4.2790, 4.7834};
static const double TDL_B_POWER_DB[TDL_B_TAPS] = {
    0.0,  -2.2, -4.0, -3.2, -9.8, -1.2,  -3.4, -5.2,  -7.6,  -3.0, -8.9, -9.0,
    -4.8, -5.7, -7.5, -1.9, -7.6, -12.2, -9.8, -11.4, -14.9, -9.2, -11.3};

/* -------------------- TDL-C -------------------- */
#define TDL_C_TAPS 24
static const double TDL_C_DELAY[TDL_C_TAPS] = {
    0.0000, 0.2099, 0.2219, 0.2329, 0.2176, 0.6366, 0.6448, 0.6560,
    0.6584, 0.7935, 0.8213, 0.9336, 1.2285, 1.3083, 2.1704, 2.7105,
    4.2589, 4.6003, 5.4902, 5.6077, 6.3065, 6.6374, 7.0427, 8.6523};
static const double TDL_C_POWER_DB[TDL_C_TAPS] = {
    -4.4,  -1.2,  -3.5,  -5.2,  -2.5,  0.0,   -2.2,  -3.9,
    -7.4,  -7.1,  -10.7, -11.1, -5.1,  -6.8,  -8.7,  -13.2,
    -13.9, -13.9, -15.8, -17.1, -16.0, -15.7, -21.6, -22.8};

/* ------------------------------------------------------------
 * TDL Channel Object
 * ------------------------------------------------------------ */
typedef struct {
  tdl_type_t type;  /**< TDL-A/B/C */
  double fs;        /**< Sampling rate [Hz] */
  double fc;        /**< Carrier frequency [Hz] */
  double speed_kmh; /**< Terminal speed [km/h] */

  int num_taps; /**< Number of taps */

  double *delay_norm; /**< Normalized delay τ_model,n (dimensionless) */
  double *power_lin;  /**< Linear power = 10^(dB/10) */

  complex float *h; /**< Rayleigh fading taps (updated each time step) */
} tdl_channel_t;

/* ------------------------------------------------------------
 * API
 * ------------------------------------------------------------ */

/**
 * @brief Create TDL channel.
 *
 * @param type       TDL_A / TDL_B / TDL_C
 * @param fs         Sampling rate [Hz]
 * @param fc         Carrier frequency [Hz]
 * @param speed_kmh  Terminal speed [km/h]
 * @return Pointer to tdl_channel_t
 */
tdl_channel_t *tdl_create(tdl_type_t type, double fs, double fc,
                          double speed_kmh);

/** @brief Free TDL channel object. */
void tdl_free(tdl_channel_t *ch);

/**
 * @brief Update fading taps h[l] using Jakes-like AR(1) model.
 *
 * Correlation:
 *   ρ = J0(2π fd Ts)
 *   σ = sqrt(1 − ρ²)
 */
void tdl_update(tdl_channel_t *ch);

/**
 * @brief Build discrete CIR from taps and per-tap delays (in seconds).
 *
 * @param ch         Channel object
 * @param delay_sec  Array delay_sec[l] = τ_scaled[l] (seconds)
 * @param out        Output CIR (length = N)
 * @param N          Number of samples in CIR
 */
void tdl_get_channel(const tdl_channel_t *ch, const double *delay_sec,
                     complex float *out, int N);

#endif /* CHANNEL_TDL_H */
