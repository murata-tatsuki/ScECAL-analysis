include ../channel_response/root.mk
.PHONY: all
all: Aggregate Plot
Aggregate: Aggregate.cc Common.hh ../calibration_parameters/EBUdecode.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) Aggregate.cc ../calibration_parameters/EBUdecode.cxx $(ROOTLIBS) -o $@
Plot: Plot.cc Common.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) Plot.cc $(ROOTLIBS) -o $@
