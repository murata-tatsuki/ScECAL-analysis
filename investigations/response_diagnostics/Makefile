CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17 -I../../analysis/impact_studies
ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
PROGRAMS = SurveyImpact AuditMIP InspectResponse InspectTopology CheckTruthResponse PlotResponse
.PHONY: all clean
all: $(PROGRAMS)
$(PROGRAMS): %: %.cc ../../analysis/impact_studies/StudyCommon.hh ../../analysis/impact_studies/EventHistograms.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
clean:
	rm -f $(PROGRAMS)
