CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17
CPPFLAGS = -I../sipm_gain/common
ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
BUILD_TMP := $(abspath ../../.scratch/sipm_gain_pool_build)
HEADERS = ../sipm_gain_selection/Common.hh ../sipm_gain/compare/CompareIO.hh $(wildcard ../sipm_gain/common/*.hh)
.PHONY: all clean
all: Pool
Pool: Pool.cc $(HEADERS) Makefile
	@set -eu; mkdir -p "$(BUILD_TMP)"; tmp=$$(mktemp "$(BUILD_TMP)/Pool.XXXXXX"); \
	trap 'rm -f -- "$$tmp"' EXIT; TMPDIR="$(BUILD_TMP)" $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(ROOTFLAGS) Pool.cc -o "$$tmp" $(ROOTLIBS) -lcrypto; \
	chmod 755 "$$tmp"; mv -f -- "$$tmp" "$@"; trap - EXIT; rmdir "$(BUILD_TMP)" 2>/dev/null || :
clean:
	rm -f Pool
