include ../../analysis/channel_response/root.mk
CXXFLAGS += -I../../analysis/impact_studies
TruthResponse: TruthResponse.cc ../../analysis/impact_studies/StudyCommon.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
Sensitivity: Sensitivity.cc ../../analysis/impact_studies/StudyCommon.hh ../../analysis/resolution/HLTailSelection.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
PrepareReference: PrepareReference.cc ../../analysis/impact_studies/StudyCommon.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
ResolutionImpact: ResolutionImpact.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
IntegratedReview: IntegratedReview.cc CurrentParameters.hh ResolutionImpact.cc ../../analysis/resolution/HLTailSelection.hh ../../analysis/impact_studies/StudyCommon.hh
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
