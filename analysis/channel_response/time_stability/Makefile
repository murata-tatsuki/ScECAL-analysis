include ../root.mk
DATA_ANALYSIS ?= /home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO
GEOMETRY = $(DATA_ANALYSIS)/src/EBUdecode.cxx
HEADERS = Common.hh ../muon_mip_response/MuonMIPCommon.hh ../muon_mip_response/TrackGeometry.hh ../muon_mip_response/CurrentCalibration.hh
BUILD_TMP = $(abspath ../../../.scratch/time_stability_build)
define build
	@set -eu; mkdir -p "$(BUILD_TMP)"; tmp=$$(mktemp "$(BUILD_TMP)/$@.XXXXXX"); \
	trap 'rm -f -- "$$tmp"' EXIT HUP INT TERM; \
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -I$(DATA_ANALYSIS)/include $< $(GEOMETRY) -o "$$tmp" $(ROOTLIBS) $(1); \
	chmod 755 "$$tmp"; mv -f -- "$$tmp" "$@"; \
	trap - EXIT HUP INT TERM; rmdir "$(BUILD_TMP)" 2>/dev/null || :
endef
all: Collect Analyze Draw
Collect: Collect.cc $(HEADERS) ../../resolution/HLTailSelection.hh $(GEOMETRY)
	$(call build,)
Analyze: Analyze.cc DiskReduce.hh StreamingStats.hh HistogramStore.hh FitReuse.hh FitPool.hh $(HEADERS) $(GEOMETRY)
	$(call build,-lz)
Draw: Draw.cc DiskReduce.hh $(HEADERS) ../adc_energy/FastPng.hh $(GEOMETRY)
	$(call build,-lASImage)
.PHONY: all
