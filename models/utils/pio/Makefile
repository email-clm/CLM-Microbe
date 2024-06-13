#
# Defined externally in Makefile.conf
#
#   INCLUDES
#   LIBS
#   MPICC
#   MPIFC
#   COPTS
#   FOPTS
#   CFLAGS
#   FFLAGS
#   AWK
#   AR

RM=/bin/rm -f
FDEPENDS=$(SRCDIR)/fdepends.awk
DEPSUF = .d
MODSUF = .mod
CPPSUF = .f90
GENF90 := ../bin

ifeq (,$(PIOARCH))
  PIOARCH=conf
endif
include Makefile.$(PIOARCH)
export PIOARCH

ifeq ($(findstring BG,$(CPPDEFS)),BG)
  SRCS_C = topology.c
else
  SRCS_C =  
endif

SRCS_FC =  pio.F90 \
            pio_kinds.F90  \
            nf_mod.F90     \
            ionf_mod.F90 \
            pio_types.F90  \
            calcdecomp.F90 \
            piolib_mod.F90 \
	    pio_mpi_utils.F90 \
            pio_nf_utils.F90 \
	    pio_utils.F90 \
            pio_msg_mod.F90 \
	    calcdisplace_mod.F90 \
	    pio_msg_callbacks.F90 \
            pio_support.F90 \
            C_interface_mod.F90

ifeq (yes,$(USE_COMPRESSION))
SRCS_FC += piovdc.F90
endif

TEMPLATES_FC = pionfatt_mod.F90.in \
	        pionfread_mod.F90.in \
                pionfwrite_mod.F90.in \
                pionfput_mod.F90.in \
                pionfget_mod.F90.in \
	        alloc_mod.F90.in \
                box_rearrange.F90.in \
                rearrange.F90.in \
	        iompi_mod.F90.in  \
	        piodarray.F90.in \
	        pio_spmd_utils.F90.in \
	        pio_msg_getput_callbacks.F90.in

TEMPSRCFC = $(TEMPLATES_FC:.in=)
SRCS_FC += $(TEMPSRCFC)
MOD_SRCS := $(filter-out C_interface_mod.F90 pio_msg_callbacks.F90 pio_msg_getput_callbacks.F90, $(SRCS_FC))

OBJS=  $(SRCS_C:.c=.o) \
       $(SRCS_FC:.F90=.o)

PP_SRC=$(SRCS_FC:.F90=.f90)


MODFILES := $(MOD_SRCS:.F90=$(MODSUF))

# File foo.DEPSUF will contain make dependency for foo.F90
DEPENDS := $(SRCS_FC:.F90=$(DEPSUF))
PERL = /usr/bin/perl

LIB= libpio.a

LEVEL=0


all: depends
	@$(MAKE) LEVEL=1 $(LIB)



depends: $(DEPENDS)
	@echo "Done updating dependencies"

ppsrc: $(PP_SRC)


#
# Suffix rules
#

.SUFFIXES:
.SUFFIXES: .o .c .F90  $(DEPSUF)


# For each file foo.F90, make the depency file foo.d

%$(DEPSUF): %.F90
	@echo 'Making dependencies for' $< '-->' $@
	@$(AWK) -f $(FDEPENDS) -v NAME=$(basename $<) -v SUF=$(suffix $<) $< > $@


SRCS_CPP=
%.o : %.F90
	$(MPIFC) -c $(FPPDEFS) $(FFLAGS) $(FOPTS) $(INCLUDES) $<

%.o : %.c
	$(MPICC) -c $(CPPDEFS) $(CFLAGS) $(COPTS) $(INCLUDES) $<

%.f90 : %.F90
	$(CPP) $< -o  $@

$(TEMPSRCFC): 
	$(PERL) $(GENF90)/genf90.pl $< > $*.F90



$(LIB): $(OBJS)
	$(RM) $@
	$(AR) $(ARFLAGS) $@ $(OBJS)


predist: predistclean $(TEMPSRCFC)

clean:
	$(RM) $(LIB) $(OBJS) $(MODFILES) $(DEPENDS) $(SRCS_CPP)

predistclean: clean
	$(RM) $(TEMPSRCFC)

install: $(LIB)  
	$(INSTALL) -d -m 755 $(INSTALL_LIB)
	$(INSTALL) -d -m 755 $(INSTALL_INC)
	$(INSTALL) -m 644 $(LIB) $(INSTALL_LIB)/$(LIB)
	$(INSTALL) -m 644 $(MODFILES) $(INSTALL_INC)


startandcount: calcdecomp.F90
	$(FC) $(FC_DEFINE)TESTCALCDECOMP $(FFLAGS) $< -o $@

#
# Automatically generated module dependencies
#

ifneq (0,$(LEVEL))
  include $(DEPENDS)
endif

alloc_mod.o : alloc_mod.F90
box_rearrange.o : box_rearrange.F90
iompi_mod.o : iompi_mod.F90
piodarray.o : piodarray.F90
pionfatt_mod.o : pionfatt_mod.F90
pionfget_mod.o : pionfget_mod.F90
pionfput_mod.o : pionfput_mod.F90
pionfread_mod.o : pionfread_mod.F90
pionfwrite_mod.o : pionfwrite_mod.F90
pio_support.o : pio_support.F90
rearrange.o : rearrange.F90

alloc_mod.F90 : alloc_mod.F90.in
box_rearrange.F90 : box_rearrange.F90.in
iompi_mod.F90 : iompi_mod.F90.in
piodarray.F90 : piodarray.F90.in
pionfatt_mod.F90 : pionfatt_mod.F90.in
pionfget_mod.F90 : pionfget_mod.F90.in
pionfput_mod.F90 : pionfput_mod.F90.in
pionfread_mod.F90 : pionfread_mod.F90.in
pionfwrite_mod.F90 : pionfwrite_mod.F90.in
rearrange.F90 : rearrange.F90.in
pio_spmd_utils.F90: pio_spmd_utils.F90.in
pio_msg_getput_callbacks.F90 : pio_msg_getput_callbacks.F90.in

