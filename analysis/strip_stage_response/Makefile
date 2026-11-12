include ../channel_response/root.mk
.PHONY: all
all: Extract Plot
Extract: Extract.cc ../ssa_shower_shape/Common.hh ../calibration_parameters/EBUdecode.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) Extract.cc ../calibration_parameters/EBUdecode.cxx $(ROOTLIBS) -o $@
Plot: Plot.cc ../beam_quantitative/PlotCommon.hh ../calibration_parameters/EBUdecode.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) Plot.cc ../calibration_parameters/EBUdecode.cxx $(ROOTLIBS) -o $@
