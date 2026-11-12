include ../root.mk
PROGRAMS = InspectMuonInputs CollectMuonMIP MergeMuonShards AnalyzeMuonMIP RedrawMuonMIP
GEOMETRY = ../../calibration_parameters/EBUdecode.cxx
BUILD_TMP = $(abspath ../../../.scratch/muon_mip_response_build)
# Replace executables atomically: running jobs keep their existing binary inode.
define build_program
	@set -eu; mkdir -p "$(BUILD_TMP)"; \
	tmp=$$(mktemp "$(BUILD_TMP)/$(notdir $@).XXXXXX"); \
	trap 'rm -f -- "$$tmp"' EXIT HUP INT TERM; \
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(1) -o "$$tmp" $(ROOTLIBS) $(2); \
	chmod 755 "$$tmp"; mv -f -- "$$tmp" "$@"; \
	trap - EXIT HUP INT TERM; rmdir "$(BUILD_TMP)" 2>/dev/null || :
endef
all: $(PROGRAMS)
InspectMuonInputs: InspectMuonInputs.cc MuonMIPCommon.hh
	$(call build_program,,)
CollectMuonMIP: CollectMuonMIP.cc MuonMIPCommon.hh TrackGeometry.hh $(GEOMETRY)
	$(call build_program,$(GEOMETRY),)
MergeMuonShards: MergeMuonShards.cc MuonMIPCommon.hh MuonMIPPlotStyle.hh ../adc_energy/FastPng.hh
	$(call build_program,,-lASImage)
AnalyzeMuonMIP: AnalyzeMuonMIP.cc CurrentCalibration.hh MuonMIPCommon.hh ExactFitCache.hh MuonMIPPlotStyle.hh ../adc_energy/FastPng.hh $(GEOMETRY) MergeMuonShards
	$(call build_program,$(GEOMETRY),-lASImage)
RedrawMuonMIP: RedrawMuonMIP.cc MuonMIPCommon.hh MuonMIPPlotStyle.hh ../adc_energy/FastPng.hh
	$(call build_program,,-lASImage)
clean:
	rm -f -- $(PROGRAMS)
.PHONY: all clean
