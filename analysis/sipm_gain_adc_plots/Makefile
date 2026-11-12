CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17
ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
BUILD_TMP := $(abspath ../../.scratch/sipm_gain_adc_plots_build)
.PHONY: all clean
all: Plot
Plot: Plot.cc FftReconstruction.hh Makefile ../sipm_gain/common/GainIO.hh ../sipm_gain/common/GainMethods.hh
	@set -eu; mkdir -p "$(BUILD_TMP)"; tmp=$$(mktemp "$(BUILD_TMP)/Plot.XXXXXX"); \
	trap 'rm -f -- "$$tmp"' EXIT; trap 'exit 143' HUP INT TERM; \
	TMPDIR="$(BUILD_TMP)" $(CXX) $(CXXFLAGS) $(ROOTFLAGS) Plot.cc -o "$$tmp" $(ROOTLIBS) -lcrypto; \
	chmod 755 "$$tmp"; mv -f -- "$$tmp" "$@"; rmdir "$(BUILD_TMP)"; trap - EXIT HUP INT TERM
clean:
	rm -f Plot
