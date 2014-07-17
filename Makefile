#
# Makefile to build ERSEM models
#

include ../../../Rules.make

DOCSRC	= shared.F90 pelagic_base.F90 benthic_base.F90 oxygen.F90 carbonate.F90 light.F90 primary_producer.F90 microzooplankton.F90 mesozooplankton.F90 bacteria.F90 calcification.F90 benthic_nutrients.F90

# Objects that derive from pelagic_base:
PELAGIC_OBJS = ${LIBFABM}(primary_producer.o) ${LIBFABM}(microzooplankton.o) ${LIBFABM}(mesozooplankton.o) ${LIBFABM}(bacteria.o) ${LIBFABM}(calcification.o)

# Objects that derive from benthic_base:
BENTHIC_OBJS = 

# All objects
OBJS = ${LIBFABM}(pelagic_base.o) ${LIBFABM}(benthic_base.o) ${LIBFABM}(oxygen.o) ${LIBFABM}(carbonate.o) ${LIBFABM}(light.o) ${LIBFABM}(benthic_nutrients.o) ${PELAGIC_OBJS} ${BENTHIC_OBJS}

all: objs

${PELAGIC_OBJS}: ${LIBFABM}(pelagic_base.o)

${BENTHIC_OBJS}: ${LIBFABM}(benthic_base.o)

${OBJS}: ${LIBFABM}(shared.o) $(FABMBASE)

${LIBFABM}(shared.o): $(FABMBASE)

objs: ${OBJS}

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/pml_ersem.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Jorn Bruggeman (BB)                             !
#-----------------------------------------------------------------------
