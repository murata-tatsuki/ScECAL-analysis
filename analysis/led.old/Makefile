TARGET = selection histo_fft selection_mip
CC = g++

CFLAGS = -O2 -Wall -Wno-deprecated
FLAG := -I.
LIBS :=
ROOTFLAG := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

.SUFFIXES: .cc .o

.cc.o :
	@echo "Now generating object files, wait a minutes..."
	$(CC) $(CFLAGS) $(FLAG) $(ROOTFLAG) -c $<

all: $(TARGET)

selection : selection.o
	@echo "Now generating target, wait a minutes..."
	$(CC) $(CFLAGS) $(LIBS) $(ROOTLIBS) $< -o $@

histo_fft : histo_fft.o
	@echo "Now generating target, wait a minutes..."
	$(CC) $(CFLAGS) $(LIBS) $(ROOTLIBS) $< -o $@

selection_mip : selection_mip.o
	@echo "Now generating target, wait a minutes..."
	$(CC) $(CFLAGS) $(LIBS) $(ROOTLIBS) $< -o $@

clean :
	rm *.o 
cleanall :
	rm *~ $(OBJS) $(TARGET)
