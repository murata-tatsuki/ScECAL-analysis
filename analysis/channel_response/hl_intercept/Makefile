include ../root.mk
.PHONY: all test clean
all: CompareHLIntercept
CompareHLIntercept: CompareHLIntercept.cc ../root.mk
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
TestComparison: tests/TestComparison.cc CompareHLIntercept.cc ../root.mk
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@
test: TestComparison
	@set -e; tmp=$$(mktemp -d); trap 'rm -rf -- "$$tmp"' EXIT; ./TestComparison "$$tmp"
clean:
	rm -f CompareHLIntercept TestComparison
