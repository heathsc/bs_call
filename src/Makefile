#==================================================================================================
# PROJECT: bscall
# FILE: Makefile
# DATE: 02/05/2017
# AUTHOR(S): Marcos Fernandez <marcos.fernandez@cnag.crg.eu>
# DESCRIPTION: Builds the BsCall
#==================================================================================================

# Definitions
ROOT_PATH=..
include Makefile.mk

PROG=bs_call

SRC=bs_call.c band_align_ss2.c 
INC=$(FOLDER_INCLUDE)/bs_call.h

BIN=$(addprefix $(FOLDER_BIN)/, $(PROG))
LIBS:=-lgsl -lgslcblas -lpthread -lm
ifeq ($(HAVE_OPENMP),1)
LIBS:=$(LIBS) -fopenmp
endif
ifeq ($(HAVE_ZLIB),1)
LIBS:=$(LIBS) -lz
endif
ifeq ($(HAVE_BZLIB),1)
LIBS:=$(LIBS) -lbz2
endif

all: setup
all: CC_FLAGS=$(GENERAL_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZATION_FLAGS)
all: $(BIN)

debug: setup
debug: CC_FLAGS=-O0 $(GENERAL_FLAGS) $(DEBUG_FLAGS)
debug: $(BIN)

setup: 
	@mkdir -p $(FOLDER_BIN) 


$(BIN): $(SRC) $(INC)
	$(CC) $(CC_FLAGS)  -o $@ $(SRC) $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS)

clean:
	@rm $(BIN)

