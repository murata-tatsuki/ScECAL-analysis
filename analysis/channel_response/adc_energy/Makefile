include ../root.mk

.PHONY: all clean
all: ChannelResponse PlotChannelResponse

ChannelResponse: ChannelResponse.cc ../root.mk
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -o $@

# Install atomically: running plot workers may still use the previous binary.
PlotChannelResponse: PlotChannelResponse.cc FastPng.hh ../root.mk
	@set -e; target="$@.tmp.$$$$"; trap 'rm -f -- "$$target"' EXIT; \
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $< $(ROOTLIBS) -lASImage -o "$$target"; \
	mv -f -- "$$target" "$@"

clean:
	rm -f ChannelResponse PlotChannelResponse
