# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.14 2010-01-11 14:27:29 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := efficiency

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all
all: lib bin

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS = $(shell $(ROOTSYS)/bin/root-config --glibs)
EXTRALIBS += $(ROOTLIBS)

include $(G4INSTALL)/config/binmake.gmk


