CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17
ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
.PHONY: all test clean
all: HLInterceptStudy SummarizeImpact ReplotEventEnergy
HLInterceptStudy SummarizeImpact TestHLIntercept ReplotEventEnergy: %: %.cc StudyCommon.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
HLInterceptStudy ReplotEventEnergy TestHLIntercept: EventHistograms.hh
test: all TestHLIntercept
	bash test.sh
	bash test_all.sh
clean:
	rm -f HLInterceptStudy SummarizeImpact TestHLIntercept ReplotEventEnergy
