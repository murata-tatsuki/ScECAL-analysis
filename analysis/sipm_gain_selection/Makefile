CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17
CPPFLAGS = -I../sipm_gain/common
ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
BUILD_TMP := $(abspath ../../.scratch/sipm_gain_selection_build)
SOURCES = Select.cc Merge.cc Plots.cc Temperature.cc
OBJECTS = $(SOURCES:.cc=.o)
HEADERS = Common.hh Selection.hh ../sipm_gain/compare/CompareIO.hh $(wildcard ../sipm_gain/common/*.hh)
.PHONY: all clean
all: Select
	@rmdir "$(BUILD_TMP)" 2>/dev/null || :
define build
	@set -eu; mkdir -p "$(BUILD_TMP)"; tmp=$$(mktemp "$(BUILD_TMP)/$(notdir $@).XXXXXX"); \
	trap 'rm -f -- "$$tmp"' EXIT; trap 'exit 143' HUP INT TERM; \
	TMPDIR="$(BUILD_TMP)" $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(ROOTFLAGS) $(1) -o "$$tmp" $(2); \
	case "$@" in *.o) chmod 644 "$$tmp";; *) chmod 755 "$$tmp";; esac; \
	mv -f -- "$$tmp" "$@"; trap - EXIT HUP INT TERM
endef
Select: $(OBJECTS)
	$(call build,$^,$(ROOTLIBS) -lcrypto)
%.o: %.cc $(HEADERS) Makefile
	$(call build,-c $<,)
clean:
	rm -f Select $(OBJECTS)
