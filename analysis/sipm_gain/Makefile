CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17
ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
CPPFLAGS = -Icommon
BUILD_TMP := $(abspath ../../.scratch/sipm_gain_build)
SOURCES = collect/Collections.cc workflow/Analyze.cc workflow/Submit.cc compare/Compare.cc common/GainMethods.cc common/Numerics.cc common/CombFit.cc multi_gaussian/Fit.cc envelope_comb/Fit.cc fft/Fit.cc autocorrelation/Fit.cc
OBJECTS = $(SOURCES:.cc=.o)
HEADERS = $(wildcard common/*.hh)
.PHONY: all clean
all: Collect Analyze
	@rmdir "$(BUILD_TMP)" 2>/dev/null || :
define build
	@set -eu; mkdir -p "$(BUILD_TMP)"; tmp=$$(mktemp "$(BUILD_TMP)/$(notdir $@).XXXXXX"); \
	trap 'rm -f -- "$$tmp"' EXIT; trap 'exit 143' HUP INT TERM; \
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(ROOTFLAGS) $(1) -o "$$tmp" $(2); \
	case "$@" in *.o) chmod 644 "$$tmp";; *) chmod 755 "$$tmp";; esac; \
	mv -f -- "$$tmp" "$@"; trap - EXIT HUP INT TERM
endef
Collect: collect/Collect.o
	$(call build,$^,$(ROOTLIBS))
Analyze: $(OBJECTS)
	$(call build,$^,$(ROOTLIBS) -lMinuit2 -lcrypto)
%.o: %.cc $(HEADERS) Makefile
	$(call build,-c $<,)
compare/Compare.o: compare/CompareIO.hh
collect/Collect.o: ../channel_response/temperature_calibration/TemperatureGeometry.hh
clean:
	rm -f Collect Analyze collect/Collect.o $(OBJECTS)
