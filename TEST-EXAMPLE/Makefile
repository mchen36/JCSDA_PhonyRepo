############################################################################################
#
# GNU Makefile 
# This is a GNU makefile to compile the F90/F03 module files of Infrared Water surface 
# modeling system. When the system is expanded/improved, users/developers may add their
# files to the OBJ_FILE list and recompile the system.
#
# CREATION HISTORY:
#       Written by:     Ming Chen, 28-Feb-2015
#                       ming.chen@noaa.gov
#  
# 
############################################################################################

include make.macros

# PATH to the shared source files
LIB_DIR = /data/jcsda/mchen/CSEM/csem_v1.0.0

INC_FILE = -I$(LIB_DIR)/inc
LIB_FILE = $(LIB_DIR)/lib/libCSEM.a 

OBJ_FILE = 

LDFLAGS_HDF5 =  $(HDF5_LIB) 
FFLAGS_HDF5  =  $(HDF5_INC) 
EXE_FILE1 = TEST_CSEM_WaterIR_SfcOptics
EXE_FILE2 = TEST_CSEM_WaterIR_BRDF
EXE_FILE3 = TEST_CSEM_LandMW_SfcOptics
EXE_FILE4 = TEST_CSEM_WaterMW_SfcOptics

EXE_FILE = $(EXE_FILE1)

all:  $(EXE_FILE)

$(EXE_FILE): $ ${OBJ_FILE} ${EXE_FILE}.o
	$(FL) $(NETCDF_LIB) $(LDFLAGS_HDF5) $(FLFLAGS) -o $(EXE_FILE) \
	${OBJ_FILE} ${EXE_FILE}.o $(LIB_FILE) 



# A clean target
clean:
	@echo "=============================================="
	@echo "Cleaning up IR Water library files"
	@echo "=============================================="
	@echo; rm -f $(OBJ_FILE)  *.o *.mod ${EXE_FILE} 

# Suffix rule definitions
# -----------------------
.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
	$(FC) -c $(FCFLAGS) $(INC_FILE) $(NETCDF_INC) $(FFLAGS_HDF5) $<
