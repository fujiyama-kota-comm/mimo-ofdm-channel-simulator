# ============================================================================
# Makefile for MIMO-OFDM Channel Simulator
# ----------------------------------------
# Builds the BER simulator for 1×1〜2×2 MIMO-OFDM under 3GPP TDL fading.
#
# Components:
#   - OFDM processing (IFFT/FFT, CP)
#   - Multi-QAM modulator/demodulator
#   - 3GPP TDL-A/B/C fading model with Jakes Doppler
#   - MIMO equalizers (ZF / MMSE)
#
# Output binary:
#   bin/mimo_ofdm_ber  (or .exe on Windows)
#
# Usage:
#   make        → build simulator
#   make run    → execute binary
#   make clean  → remove build files
# ============================================================================

CC      = gcc
CFLAGS  = -O2 -Wall -std=c99 -Iinclude
LDFLAGS = -lm

BIN_DIR = bin

# ----------------------------------------------------------------------------
# Source files
# ----------------------------------------------------------------------------
SRC_COMMON = \
    src/modulator.c \
    src/ofdm.c \
    src/channel_tdl.c \
    src/mimo.c

OBJ_COMMON = $(SRC_COMMON:.c=.o)

# ----------------------------------------------------------------------------
# Main entry (MIMO-OFDM BER simulator)
# ----------------------------------------------------------------------------
MIMO_SRC = mains/mimo_ofdm_ber.c
MIMO_OBJ = $(MIMO_SRC:.c=.o)

ifeq ($(OS),Windows_NT)
    TARGET = $(BIN_DIR)/mimo_ofdm_ber.exe
else
    TARGET = $(BIN_DIR)/mimo_ofdm_ber
endif

# ============================================================================
# Build rules
# ============================================================================
all: $(TARGET)

$(BIN_DIR):
	@if [ ! -d "$(BIN_DIR)" ]; then \
		mkdir -p $(BIN_DIR) 2>/dev/null || mkdir $(BIN_DIR); \
	fi

$(TARGET): $(BIN_DIR) $(OBJ_COMMON) $(MIMO_OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ_COMMON) $(MIMO_OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# ============================================================================
# Run shortcut
# ============================================================================
run:
	./$(TARGET)

# ============================================================================
# Clean
# ============================================================================
clean:
	@echo "Cleaning object files..."
	rm -f $(OBJ_COMMON) $(MIMO_OBJ)

	@echo "Cleaning binary..."
	@if [ -f "$(TARGET)" ]; then rm -f "$(TARGET)"; fi

	@if [ -d "$(BIN_DIR)" ] && [ ! "$$(ls -A $(BIN_DIR))" ]; then \
		echo "Removing empty bin directory"; \
		rmdir $(BIN_DIR); \
	fi

.PHONY: all clean run
