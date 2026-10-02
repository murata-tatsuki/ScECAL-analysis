CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17
ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
.PHONY: all check clean
all: CollectMuon
CollectMuon: CollectMuon.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@.tmp
	mv $@.tmp $@
check: CollectMuon
	python3.12 test_mip_closure.py
clean:
	rm -f CollectMuon CollectMuon.tmp
