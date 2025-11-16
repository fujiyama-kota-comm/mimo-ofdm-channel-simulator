"""
plot_mimo_ofdm_ber.py
---------------------
Plot BER curves for 2×2 MIMO-OFDM over a 3GPP TDL fading channel.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# =============================================================================
#  Setup
# =============================================================================

CSV_PATH = "results/mimo_ofdm_ber_data.csv"
OUTPUT_DIR = "images"
os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14

# =============================================================================
#  Load CSV
# =============================================================================
if not os.path.exists(CSV_PATH):
    print(f"[ERROR] CSV file not found: {CSV_PATH}")
    sys.exit(1)

df = pd.read_csv(CSV_PATH)

SNR = df["SNR_dB"]

mod_cols = [c for c in df.columns if c.startswith("BER_") and c != "SNR_dB"]

preferred_order = ["BER_BPSK", "BER_QPSK", "BER_16QAM", "BER_64QAM", "BER_256QAM"]
mod_cols = [c for c in preferred_order if c in mod_cols]

labels = {
    "BER_BPSK": "BPSK (1 bit/sym)",
    "BER_QPSK": "QPSK (2 bit/sym)",
    "BER_16QAM": "16-QAM (4 bit/sym)",
    "BER_64QAM": "64-QAM (6 bit/sym)",
    "BER_256QAM": "256-QAM (8 bit/sym)",
}

markers = ["o", "s", "^", "v", "D"]

colors = [
    "#1f77b4",  # Blue
    "#2ca02c",  # Green
    "#9467bd",  # Purple
    "#ff7f0e",  # Orange
    "#8c564b",  # Brown
]

# =============================================================================
#  Plot
# =============================================================================

plt.figure(figsize=(8, 6))

for col, marker, color in zip(mod_cols, markers, colors):
    plt.semilogy(
        SNR,
        df[col],
        marker=marker,
        markersize=8,
        markerfacecolor="none",
        linewidth=2.5,
        color=color,
        label=labels.get(col, col),
    )

plt.ylim(1e-5, 1)
plt.xlabel("SNR [dB]", fontsize=18)
plt.ylabel("Bit Error Rate (BER)", fontsize=18)

plt.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.6)
plt.legend(fontsize=13, loc="upper right", frameon=True, edgecolor="black")
plt.tight_layout()

out_png = os.path.join(OUTPUT_DIR, "mimo_ofdm_ber_graph.png")
out_svg = os.path.join(OUTPUT_DIR, "mimo_ofdm_ber_graph.svg")

plt.savefig(out_png, dpi=300)
plt.savefig(out_svg, dpi=300)

print(f"[OK] Saved: {out_png}")
print(f"[OK] Saved: {out_svg}")

plt.show()
