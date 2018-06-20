#============================================================================
# PROJECT: bscall
# FILE: Makefile
# DATE: 02/05/2017
# AUTHOR(S): Marcos Fernandez Callejo <mfernandez@cnag.crg.eu>
# DESCRIPTION: Top level makefile
#============================================================================

# Definitions
all: release

release: 	
	$(MAKE) --directory=GEMTools
	$(MAKE) --directory=src

debug:
	$(MAKE) --directory=GEMTools debug
	$(MAKE) --directory=src debug

static:
	$(MAKE) --directory=GEMTools static
	$(MAKE) --directory=src static

clean:
	$(MAKE) --directory=GEMTools clean
	$(MAKE) --directory=src clean

