include ../root.mk
PYTHON ?= python3.12
.PHONY: all check clean
all: CollectTemperature
CollectTemperature: CollectTemperature.cc TemperatureGeometry.hh ../calibration_residual/RootInput.hh ../calibration_residual/ResidualStats.hh ../root.mk
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@.tmp
	mv $@.tmp $@
check: CollectTemperature
	$(PYTHON) test_temperature.py
clean:
	rm -f CollectTemperature CollectTemperature.tmp
