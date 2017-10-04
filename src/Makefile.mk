#==================================================================================================
# PROJECT: bs_call
# FILE: Makefile.mk
# DATE: 02/05/2017
# AUTHOR(S): Simon Heath (simon.heath@gmail.com) and Marcos Fernandez (marcos.fernandez@cnag.crg.eu)
# DESCRIPTION: Makefile definitions' file
#==================================================================================================

#Include
include ../Gsl.mk


# Utilities
CC=/usr/bin/gcc
AR=ar

# Folders
FOLDER_BIN=../bin
FOLDER_INCLUDE=../include

# Gem tools library
GEM_TOOLS = 

# Flags
ARCH_FLAGS = -D__LINUX__
HAVE_ZLIB = 1
HAVE_BZLIB = 1
HAVE_OPENMP = 1
GEMTOOLS_INC = -I../GEMTools/include -I../GEMTools/resources/include
GEMTOOLS_LIBS = -L../GEMTools/lib -lgemtools

GENERAL_FLAGS=-fPIC -Wall -msse4 -std=gnu99
ifeq ($(HAVE_ZLIB),1)
GENERAL_FLAGS:=$(GENERAL_FLAGS) -DHAVE_ZLIB
endif
ifeq ($(HAVE_BZLIB),1)
GENERAL_FLAGS:=$(GENERAL_FLAGS) -DHAVE_BZLIB
endif
ifeq ($(HAVE_OPENMP),1)
GENERAL_FLAGS:=$(GENERAL_FLAGS) -DHAVE_OPENMP
endif

OPTIMIZATION_FLAGS=-O3 -g # -fomit-frame-pointer -ftree-vectorize
ARCH_FLAGS_OPTIMIZATION_FLAGS= # -msse3 -mssse3 -msse4.2

INCLUDE_FLAGS=-I$(FOLDER_INCLUDE) $(GEMTOOLS_INC) $(GSL_INC)
LIB_PATH_FLAGS=$(GEMTOOLS_LIBS) $(GSL_LIB)

SUPPRESS_CHECKS=-DNDEBUG -DGT_NO_CONSISTENCY_CHECKS
DEBUG_FLAGS=-g -DGT_DEBUG

PLATFORM=$(shell uname)
