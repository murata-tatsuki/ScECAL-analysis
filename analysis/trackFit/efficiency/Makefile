TARGET = channel_efficiecy
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

channel_efficiecy : channel_efficiecy.o
	@echo "Now generating target, wait a minutes..."
	$(CC) $(CFLAGS) $(LIBS) $(ROOTLIBS) $< -o $@

clean :
	rm *.o 
cleanall :
	rm *~ $(OBJS) $(TARGET)
