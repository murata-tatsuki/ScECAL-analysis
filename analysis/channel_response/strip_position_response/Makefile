include ../root.mk
DATA_ANALYSIS ?= /home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO
GEOMETRY = $(DATA_ANALYSIS)/src/EBUdecode.cxx
CXXFLAGS += -I$(DATA_ANALYSIS)/include
PROGRAMS = Collect Aggregate Fit Draw
HEADERS = Common.hh ../muon_mip_response/MuonMIPCommon.hh ../muon_mip_response/CurrentCalibration.hh ../muon_mip_response/TrackGeometry.hh
BUILD_TMP = $(abspath ../../../.scratch/strip_position_response_build)
all: $(PROGRAMS)
define build
	@set -eu; mkdir -p "$(BUILD_TMP)"; \
	target=$$(mktemp "$(BUILD_TMP)/$@.XXXXXX"); \
	trap 'rm -f -- "$$target"' EXIT HUP INT TERM; \
	TMPDIR="$(BUILD_TMP)" $(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(1) -o "$$target" $(ROOTLIBS) $(2); \
	chmod 755 "$$target"; mv -f -- "$$target" "$@"; \
	trap - EXIT HUP INT TERM; rmdir "$(BUILD_TMP)" 2>/dev/null || :
endef
Collect: Collect.cc $(HEADERS) $(GEOMETRY)
	$(call build,$(GEOMETRY),)
Aggregate: Aggregate.cc $(HEADERS)
	$(call build,,)
Fit: Fit.cc $(HEADERS) ../../mip_calibration/FitModel.hh ../../mip_calibration/Common.hh
	$(call build,,-lMinuit2)
Draw: Draw.cc $(HEADERS) ../adc_energy/FastPng.hh
	$(call build,,-lASImage)
clean:
	rm -f -- $(PROGRAMS)
.PHONY: all clean
