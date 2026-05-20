# Top-level Makefile for bGFA
# - Ignore case/ sources
# - Build KSW object and WFA static libs
# - Place intermediates into ./build
# - Produce executable ./bgfatools at repo root

# Toolchain
CXX      ?= g++
CC       ?= gcc
UNAME_S  := $(shell uname -s)

# Paths
SRC_DIR   := src
BUILD_DIR := build
WFA_DIR   := $(SRC_DIR)/wfa
WFA_LIB   := $(WFA_DIR)/lib/libwfa.a
WFA_LIBCPP:= $(WFA_DIR)/lib/libwfacpp.a

# Sources
MAIN_SRC  := $(SRC_DIR)/main.cpp
KSW_SRC   := $(SRC_DIR)/ksw/ksw2_extz2_sse.c
KSW_OBJ   := $(BUILD_DIR)/ksw2_extz2_sse.o

# Target
TARGET := bgfatools

# Flags
CXXFLAGS ?= -std=c++17 -O3 -I./src -I./src/wfa -L./src/wfa/lib
CFLAGS   ?= -O2
LDFLAGS  ?=
LDLIBS   ?= -lm -lz -lwfacpp -lwfa

# Platform-specific
ifeq ($(UNAME_S),Linux)
  LDLIBS += -lrt
endif

.PHONY: all ksw wfa clean distclean debug test

all: $(TARGET)

$(TARGET): $(KSW_OBJ) wfa $(MAIN_SRC)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(MAIN_SRC) $(KSW_OBJ) $(WFA_LIBCPP) $(WFA_LIB) -o $(TARGET) $(LDFLAGS) $(LDLIBS)
	@echo "Built $(TARGET)"

# Build KSW SSE object into ./build
$(KSW_OBJ): $(KSW_SRC) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Build WFA static libraries via its own Makefile
wfa:
	$(MAKE) -C $(WFA_DIR) 

# Convenience debug build (adds -g, disables -O3)
debug: CXXFLAGS := $(filter-out -O3,$(CXXFLAGS)) -O0 -g -DINFO_DEBUG
debug: clean all

prof: CXXFLAGS := $(filter-out -O3,$(CXXFLAGS)) -O3 -g -pg
prof: clean all

# Minimal test: build then show help
test: all
	@./$(TARGET) --help 2>/dev/null | head -n 20 || true

clean:
	@rm -rf $(BUILD_DIR)
	@$(MAKE) -C $(WFA_DIR) clean || true
	@echo "Cleaned intermediates"
# 
# distclean: clean
	@rm -f $(TARGET)
	@echo "Removed $(TARGET)"

