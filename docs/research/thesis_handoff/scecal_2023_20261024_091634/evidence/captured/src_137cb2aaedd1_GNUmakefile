#############################################
#                                           #
############  makefile  #####################
#                                           #
#############################################

##set up control#######
TOP:= $(shell pwd)
DIR_OBJ:= $(TOP)/obj
DIR_BIN:= $(TOP)/bin
DIR_SRC:= $(TOP)/src
DIR_INCLUDE:= $(TOP)/include


SRC = $(wildcard $(DIR_SRC)/*.cxx) 	
OBJ = $(patsubst $(DIR_SRC)/%.cxx,$(DIR_OBJ)/%.o,$(SRC)) 
#ALL:
#	@echo $(SRC)
#	@echo $(OBJ)

##ROOT#######################################
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMinuit
ROOTGLIBS     = $(shell root-config --glibs)

CPPLIBS = $(ROOTLIBS) $(ROOTGLIBS) -ltbb
DIR_INCLUDE+=$(ROOTCFLAGS)

# #set up compilers#
CXX = g++
CPPFLAGS = -g -Wall -I$(DIR_INCLUDE)

####### Make Execuatbles#######
all: decode pedestal pedestal_pedestal calib track_fit track_fit_ track_fit_3d mip hl_ratio ssa calib_ssa simpleCalib

decode: decode.cxx $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

pedestal: pedestal.cxx $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

pedestal_pedestal: pedestal_pedestal.cxx $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

track_fit: track_fit.cxx $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

track_fit_: track_fit_.cxx $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

track_fit_3d: track_fit_3d.cxx $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

mip: mip.cxx $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

calib: calib.cxx $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

hl_ratio: hl_ratio.cxx $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

ssa: ssa.cxx $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

calib_ssa: calib_ssa.cxx $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

simpleCalib: simpleCalib.cxx  $(OBJ) #$(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

#################################
$(DIR_OBJ)/%.o:$(DIR_SRC)/%.cxx  
	$(CXX) $(CPPFLAGS)  -c $(DIR_SRC)/$(notdir $<)  -o $(DIR_OBJ)/$(notdir $@)

#	$(AR) r $(DIR_SRC)/$(notdir $@) $(DIR_OBJ)/$(notdir $@)
#	$(RM) $(DIR_OBJ)/$(notdir $@)
#@echo 

.PHONY:clean
clean: 
	rm -f $(DIR_OBJ)/*.o rm -f $(DIR_BIN)/*
