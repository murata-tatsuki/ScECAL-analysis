TARGET = between_files_energy_resolution_simulation between_files_energy_resolution between_files_cog between_files_energy_resolution_generation SingleEnergyAnalysis MultiEnergyAnalysis
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

between_files_energy_resolution_simulation : between_files_energy_resolution_simulation.o
	@echo "Now generating target, wait a minutes..."
	$(CC) $(CFLAGS) $(LIBS) $(ROOTLIBS) $< -o $@

between_files_energy_resolution : between_files_energy_resolution.o
	@echo "Now generating target, wait a minutes..."
	$(CC) $(CFLAGS) $(LIBS) $(ROOTLIBS) $< -o $@

between_files_cog : between_files_cog.o
	@echo "Now generating target, wait a minutes..."
	$(CC) $(CFLAGS) $(LIBS) $(ROOTLIBS) $< -o $@

between_files_energy_resolution_generation : between_files_energy_resolution_generation.o
	@echo "Now generating target, wait a minutes..."
	$(CC) $(CFLAGS) $(LIBS) $(ROOTLIBS) $< -o $@

SingleEnergyAnalysis : SingleEnergyAnalysis.o
	@echo "Now generating target, wait a minutes..."
	$(CC) $(CFLAGS) $(LIBS) $(ROOTLIBS) $< -o $@

MultiEnergyAnalysis : MultiEnergyAnalysis.o
	@echo "Now generating target, wait a minutes..."
	$(CC) $(CFLAGS) $(LIBS) $(ROOTLIBS) $< -o $@

clean :
	rm *.o 
cleanall :
	rm *~ $(OBJS) $(TARGET)
