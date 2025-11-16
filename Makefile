CC      = gcc
CFLAGS  = -O2 -Wall -std=c99 -Iinclude
LDFLAGS = -lm

BIN_DIR = bin

# ============================================================
# Source files
# ============================================================
SRC_COMMON = \
    src/modulator.c \
    src/ofdm.c \
    src/channel_tdl.c

OBJ_COMMON = $(SRC_COMMON:.c=.o)

# Main simulation
BER_SRC = mains/ofdm_ber.c
BER_OBJ = $(BER_SRC:.c=.o)

# ============================================================
# Target binary (OS dependent)
# ============================================================
ifeq ($(OS),Windows_NT)
    TARGET = $(BIN_DIR)/ofdm_ber.exe
else
    TARGET = $(BIN_DIR)/ofdm_ber
endif

# ============================================================
# Build rules
# ============================================================
all: $(TARGET)

$(BIN_DIR):
	@if [ ! -d "$(BIN_DIR)" ]; then \
		mkdir -p $(BIN_DIR) 2>/dev/null || mkdir $(BIN_DIR); \
	fi

$(TARGET): $(BIN_DIR) $(OBJ_COMMON) $(BER_OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ_COMMON) $(BER_OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# ============================================================
# Run shortcut
# ============================================================
run:
	./$(TARGET)

# ============================================================
# Clean
# ============================================================
clean:
	@echo "Cleaning object files..."
	rm -f $(OBJ_COMMON) $(BER_OBJ)

	@echo "Cleaning binaries..."
	@if [ -f "$(TARGET)" ]; then rm -f "$(TARGET)"; fi

	@if [ -d "$(BIN_DIR)" ] && [ ! "$$(ls -A $(BIN_DIR))" ]; then \
		echo "Removing empty bin directory"; \
		rmdir $(BIN_DIR); \
	fi

.PHONY: all clean run
