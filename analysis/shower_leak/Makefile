include ../channel_response/root.mk
.PHONY: all
all: Aggregate Plot
Aggregate Plot: %: %.cc Common.hh ../resolution/HLTailSelection.hh ../impact_studies/StudyCommon.hh ../calibration_parameters/EBUdecode.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< ../calibration_parameters/EBUdecode.cxx $(ROOTLIBS) -o $@

Plot: IntervalPlots.hh TrackSelection.hh PlotStyle.hh DisplayPlots.hh

