include ../root.mk
CXXFLAGS += -I../../impact_studies
GEOMETRY = ../../calibration_parameters/EBUdecode.cxx
.PHONY: all clean
all: HGLGTailStudy PlotHGLGTailLayers
HGLGTailStudy PlotHGLGTailLayers: %: %.cc HGLGTailLayerMaps.hh ../../impact_studies/StudyCommon.hh $(GEOMETRY) ../../calibration_parameters/EBUdecode.h ../root.mk
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(GEOMETRY) $(ROOTLIBS) -o $@
clean:
	rm -f HGLGTailStudy PlotHGLGTailLayers
