import pandas as pd
import matplotlib.pyplot as plt
import os

# =============================================================================
#  Plot BER Curves for OFDM under TDL Channel
# =============================================================================

os.makedirs("images", exist_ok=True)

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14

# Load CSV
df = pd.read_csv("results/ofdm_ber_data.csv")

SNR = df["SNR_dB"]

ber_bpsk = df["BER_1QAM"]
ber_qpsk = df["BER_2QAM"]
ber_16qam = df["BER_4QAM"]
ber_64qam = df["BER_6QAM"]
ber_256qam = df["BER_8QAM"]

# =============================================================================
#  Plot setup
# =============================================================================
plt.figure(figsize=(8, 6))

plt.semilogy(
    SNR,
    ber_bpsk,
    marker="o",
    markersize=8,
    markerfacecolor="none",
    linewidth=2.5,
    label="BPSK (1 bit/sym)",
)

plt.semilogy(
    SNR,
    ber_qpsk,
    marker="s",
    markersize=8,
    markerfacecolor="none",
    linewidth=2.5,
    label="QPSK (2 bit/sym)",
)

plt.semilogy(
    SNR,
    ber_16qam,
    marker="^",
    markersize=8,
    markerfacecolor="none",
    linewidth=2.5,
    label="16-QAM (4 bit/sym)",
)

plt.semilogy(
    SNR,
    ber_64qam,
    marker="v",
    markersize=8,
    markerfacecolor="none",
    linewidth=2.5,
    label="64-QAM (6 bit/sym)",
)

plt.semilogy(
    SNR,
    ber_256qam,
    marker="D",
    markersize=8,
    markerfacecolor="none",
    linewidth=2.5,
    label="256-QAM (8 bit/sym)",
)

# Y-axis limits
plt.ylim(1e-5, 1)

# Labels
plt.xlabel("SNR [dB]", fontsize=18)
plt.ylabel("Bit Error Rate (BER)", fontsize=18)

plt.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.6)
plt.legend(fontsize=13, loc="upper right", frameon=True, edgecolor="black")

plt.tight_layout()

# Save output
plt.savefig("images/ofdm_ber_graph.png", dpi=300)
plt.savefig("images/ofdm_ber_graph.svg", dpi=300)

plt.show()
