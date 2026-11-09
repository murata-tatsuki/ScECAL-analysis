TARGET = MIP_calibration
CC = g++

CFLAGS = -O2 -Wall -Wno-deprecated
FLAG := -I.
LIBS :=
ROOTFLAG := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

.PHONY: all clean cleanall
all: $(TARGET)

MIP_calibration.o: MIP_calibration.cc FitMIP.C EBUdecode.h EBUdecode.cxx
	$(CC) $(CFLAGS) $(FLAG) $(ROOTFLAG) -c $< -o $@

MIP_calibration: MIP_calibration.o
	$(CC) $(CFLAGS) $(LIBS) $< $(ROOTLIBS) -o $@

clean:
	rm -f MIP_calibration.o

cleanall: clean
	rm -f $(TARGET)
