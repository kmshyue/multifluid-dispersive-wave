#
# Makefile for multiphase code
#
FFLAGS = -O  
LFLAGS =               
FC     = gfortran     
#
OBJECTS = \
driver.o \
setprob.o \
qinit.o \
../../lib/setaux2.o \
../../lib/setaux2_grav.o \
../../lib/state2_rotate3.o \
../../lib/prmtoc.o \
../../lib/ctomfd.o \
../../lib/rpn2.o \
../../lib/rpt2.o \
../../lib/hllc.o \
../../lib/roe.o \
../../lib/src2.o \
../../lib/out2.o \
../../lib/out1.o \
../../lib/bc2.o \
../../lib/estdt2.o \
../../lib/limiter.o \
../../lib/claw2ez.o \
../../lib/claw2.o \
../../lib/claw2_maloc.o \
../../lib/step2.o \
../../lib/step2ds.o \
../../lib/dimsp2.o \
../../lib/flux2.o \
../../lib/an2step.o \
../../lib/a4step2.o \
../../lib/b4step2.o \
../../lib/mapc2p.o \
../../lib/cellave2.o \
../../lib/fss.o \
../../lib/copyaux2.o \
../../lib/copyq2.o \
../../lib/restart2.o \
../../lib/outgrd.o \
../../lib/rhoe_prho_ideal.o \
../../lib/p_rhoe_ideal.o \
../../lib/rhoe_prho_sg.o \
../../lib/p_rhoe_sg.o \
../../lib/philim.o \
../../lib/zero1d_aux.o \
../../lib/zeroin.o

SOURCES = \

xclaw: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) -o xclaw 

clear: $(OBJECTS)
	rm $(OBJECTS) xclaw

program : $(SOURCES)
	cat $(SOURCES) > claw2program.f

### DO NOT remove this line - make depends on it ###
