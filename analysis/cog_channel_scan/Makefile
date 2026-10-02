TARGET = SingleEnergyAnalysis MultiEnergyAnalysis
CXX = g++
CXXFLAGS = -O2 -Wall -Wno-deprecated
ROOTFLAG := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

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
