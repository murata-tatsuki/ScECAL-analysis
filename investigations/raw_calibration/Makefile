CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17 -I../../analysis/impact_studies
ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
all: InspectRawCalibration
InspectRawCalibration: InspectRawCalibration.cc ../../analysis/impact_studies/StudyCommon.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
all: InspectTails
InspectTails: InspectTails.cc ../../analysis/impact_studies/StudyCommon.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
all: SummarizeRaw
SummarizeRaw: SummarizeRaw.cc ../../analysis/impact_studies/StudyCommon.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
