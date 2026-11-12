include ../channel_response/root.mk
DATA_ANALYSIS ?= /home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO
GEOMETRY = $(DATA_ANALYSIS)/src/EBUdecode.cxx
CXXFLAGS += -I$(DATA_ANALYSIS)/include
PROGRAMS = Collect Fit Export Draw Diagnose
HEADERS = Common.hh ../channel_response/muon_mip_response/MuonMIPCommon.hh ../channel_response/muon_mip_response/CurrentCalibration.hh ../channel_response/muon_mip_response/TrackGeometry.hh
BUILD_TMP = $(abspath ../../.scratch/mip_calibration_build)
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
Fit: Fit.cc FitModel.hh $(HEADERS)
	$(call build,,-lMinuit2)
Diagnose: Diagnose.cc FitModel.hh $(HEADERS) ../channel_response/adc_energy/FastPng.hh
	$(call build,,-lMinuit2 -lASImage)
Export: Export.cc $(HEADERS)
	$(call build,,)
Draw: Draw.cc $(HEADERS) ../channel_response/adc_energy/FastPng.hh
	$(call build,,-lASImage)
clean:
	rm -f -- $(PROGRAMS)
.PHONY: all clean
