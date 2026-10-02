include ../root.mk
.PHONY: all clean check
all: CalibrationResidual
CalibrationResidual: CalibrationResidual.cc RootInput.hh CalibrationAudit.hh ResidualStats.hh ../adc_energy/FastPng.hh ../root.mk
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -lASImage -o $@.tmp
	mv $@.tmp $@
check: CalibrationResidual
	python3.12 test_residuals.py ./CalibrationResidual
clean:
	rm -f CalibrationResidual CalibrationResidual.tmp
