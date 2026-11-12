include ../channel_response/root.mk
.PHONY: all
all: Extract Plot ExtractThreshold PlotThreshold PlotThresholdScan
Extract: Extract.cc ../ssa_shower_shape/Common.hh ../calibration_parameters/EBUdecode.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) Extract.cc ../calibration_parameters/EBUdecode.cxx $(ROOTLIBS) -o $@
Plot: Plot.cc PlotCommon.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) Plot.cc $(ROOTLIBS) -o $@

ExtractThreshold: ExtractThreshold.cc ../ssa_shower_shape/Common.hh ../calibration_parameters/EBUdecode.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) ExtractThreshold.cc ../calibration_parameters/EBUdecode.cxx $(ROOTLIBS) -o $@
PlotThreshold: PlotThreshold.cc PlotCommon.hh ThresholdPlotCommon.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) PlotThreshold.cc $(ROOTLIBS) -o $@
PlotThresholdScan: PlotThresholdScan.cc PlotCommon.hh ThresholdPlotCommon.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) PlotThresholdScan.cc $(ROOTLIBS) -o $@
