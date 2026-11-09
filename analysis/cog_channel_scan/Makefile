TARGET = SingleEnergyAnalysis MultiEnergyAnalysis
.DEFAULT_GOAL := all
CXX = g++
CXXFLAGS = -O2 -Wall -Wno-deprecated
ROOTFLAG := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

SingleEnergyAnalysis.o MultiEnergyAnalysis.o: ../resolution/HLTailSelection.hh

.PHONY: all clean cleanall
all: $(TARGET)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAG) -c $< -o $@

$(TARGET): %: %.o
	$(CXX) $(CXXFLAGS) $< $(ROOTLIBS) -o $@

clean:
	rm -f $(TARGET:=.o)

cleanall: clean
	rm -f $(TARGET)
