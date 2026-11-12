include ../root.mk
CXXFLAGS += -I../../impact_studies
CXX_RUNTIME_DIR := $(dir $(shell env -u LIBRARY_PATH $(CXX) -print-file-name=libstdc++.so.6))
ROOTLIBS += -L$(CXX_RUNTIME_DIR) -Wl,-rpath,$(CXX_RUNTIME_DIR)
GEOMETRY = ../../calibration_parameters/EBUdecode.cxx
.PHONY: all clean
all: HGLGTailStudy PlotHGLGTailLayers
HGLGTailStudy PlotHGLGTailLayers: %: %.cc HGLGTailLayerMaps.hh HGLGTailCalibration.hh ../../impact_studies/StudyCommon.hh $(GEOMETRY) ../../calibration_parameters/EBUdecode.h ../root.mk Makefile
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(GEOMETRY) $(ROOTLIBS) -o $@
clean:
	rm -f HGLGTailStudy PlotHGLGTailLayers
