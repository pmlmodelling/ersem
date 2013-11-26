#
# Makefile to build the passive tracer model
#

include ../../../Rules.make

DOCSRC	= base.F90 nutrients.F90 dom.F90 pom.F90 gas_dynamics.F90 vphyt.F90 microzoo.F90 mesozoo.F90 bacteria.F90

OBJS    = ${LIBFABM}(base.o) ${LIBFABM}(nutrients.o) ${LIBFABM}(dom.o) ${LIBFABM}(pom.o) ${LIBFABM}(gas_dynamics.o) ${LIBFABM}(vphyt.o) ${LIBFABM}(microzoo.o) ${LIBFABM}(mesozoo.o) ${LIBFABM}(bacteria.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/pml_ersem.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Jorn Bruggeman (BB)                             !
#-----------------------------------------------------------------------
