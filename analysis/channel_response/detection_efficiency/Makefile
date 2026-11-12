include ../root.mk
PROGRAMS = Collect Analyze Draw
GEOMETRY = ../../calibration_parameters/EBUdecode.cxx
HEADERS = Common.hh ../muon_mip_response/TrackGeometry.hh ../muon_mip_response/MuonMIPCommon.hh
BUILD_TMP = $(abspath ../../../.scratch/detection_efficiency_build)
all: $(PROGRAMS)
define build
	@set -eu; mkdir -p "$(BUILD_TMP)"; tmp=$$(mktemp "$(BUILD_TMP)/$@.XXXXXX"); trap 'rm -f -- "$$tmp"' EXIT HUP INT TERM; TMPDIR="$(BUILD_TMP)" $(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(GEOMETRY) -o "$$tmp" $(ROOTLIBS) $(1); chmod 755 "$$tmp"; mv -f -- "$$tmp" "$@"; trap - EXIT HUP INT TERM; rmdir "$(BUILD_TMP)" 2>/dev/null || :
endef
Collect: Collect.cc $(HEADERS) $(GEOMETRY)
	$(call build,)
Analyze: Analyze.cc $(HEADERS) $(GEOMETRY)
	$(call build,)
Draw: Draw.cc $(HEADERS) $(GEOMETRY) ../adc_energy/FastPng.hh
	$(call build,-lASImage -lHistPainter)
clean:
	rm -f -- $(PROGRAMS)
.PHONY: all clean
