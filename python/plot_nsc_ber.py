import pandas as pd
import matplotlib.pyplot as plt
import os

# =============================================================================
#  Plot BER Curves for NSC (Non-Systematic Convolutional Code)
#  ---------------------------------------------------------------------------
#  - Input : results/nsc_ber_data.csv
#        Columns = EbN0_dB, BER_soft, BER_hard, BER_bpsk
#  - Output:
#        images/nsc_ber_graph.png
#        images/nsc_ber_graph.svg
#
#  Description:
#        Visualizes the BER performance of the NSC (rate-1/2) convolutional code
#        using soft-/hard-decision Viterbi decoding and compares it with
#        uncoded BPSK theoretical performance.
# =============================================================================

# Create output directory (if missing)
os.makedirs("images", exist_ok=True)

# =============================================================================
#  Matplotlib font settings (publication quality)
# =============================================================================
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14

# =============================================================================
#  Load CSV data
# =============================================================================
df = pd.read_csv("results/nsc_ber_data.csv")

EbN0 = df["EbN0_dB"]
ber_soft = df["BER_soft"]
ber_hard = df["BER_hard"]
ber_bpsk = df["BER_bpsk"]

# =============================================================================
#  Create BER figure
# =============================================================================
plt.figure(figsize=(7.5, 6))

# --- Soft-decision Viterbi (green)
plt.semilogy(
    EbN0,
    ber_soft,
    marker="o",
    markersize=8,
    markerfacecolor="none",
    markeredgewidth=1.8,
    linewidth=2.5,
    color="g",
    label="Soft-decision Viterbi",
)

# --- Hard-decision Viterbi (blue)
plt.semilogy(
    EbN0,
    ber_hard,
    marker="s",
    markersize=8,
    markerfacecolor="none",
    markeredgewidth=1.8,
    linewidth=2.5,
    color="b",
    label="Hard-decision Viterbi",
)

# --- Uncoded BPSK (theory, red)
plt.semilogy(
    EbN0,
    ber_bpsk,
    linewidth=2.5,
    color="r",
    label="Uncoded BPSK (theory)",
)

# Y-axis
plt.ylim(1e-5, 1)

# Axis labels
plt.xlabel("Eb/N0 [dB]", fontsize=18)
plt.ylabel("Bit Error Rate (BER)", fontsize=18)

# Grid for log-scale
plt.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.6)

# Legend
plt.legend(fontsize=14, loc="upper right", frameon=True, edgecolor="black")

# =============================================================================
#  Annotation box (code parameters)
# =============================================================================
text = (
    "NSC convolutional code\n"
    "Rate = 1/2\n"
    "Constraint length = 3\n"
    "States = 4 (2-bit register)"
)

plt.annotate(
    text,
    xy=(0.03, 0.03),
    xycoords="axes fraction",
    ha="left",
    va="bottom",
    fontsize=12,
    bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.4"),
)

plt.tight_layout()

# =============================================================================
#  Save figure (PNG + SVG)
# =============================================================================
plt.savefig("images/nsc_ber_graph.png", dpi=300, bbox_inches="tight")
plt.savefig("images/nsc_ber_graph.svg", dpi=300, bbox_inches="tight")

plt.show()
