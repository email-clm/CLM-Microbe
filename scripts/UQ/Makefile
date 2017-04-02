#CADES
#F90 = /software/dev_tools/swtree/cs400_centos7.2_pe2016-08/openmpi/1.10.3/centos7.2_gnu5.3.0/bin/mpif90
#CC = /software/dev_tools/swtree/cs400_centos7.2_pe2016-08/openmpi/1.10.3/centos7.2_gnu5.3.0/bin/mpicc
#BASE_NETCDF=/software/user_tools/current/cades-ccsi/netcdf4/4.4.0/openmpi-1.10-gcc-5.3

#OIC
#F90=/projects/cesm/devtools/mpich-3.0.4-gcc4.8.1/bin/mpif90
#CC=/projects/cesm/devtools/mpich-3.0.4-gcc4.8.1/bin/mpicc
#BASE_NETCDF=/home/zdr/netcdf-4.0

# Define environment
ENVIRON = -DERGO
ENVIRON = -DHAVE_MPI

#Compiler options
FFLAGS = -O3 -fno-range-check 
#FFLAGS = -O3 
CFLAGS = -O3
MODEL = CLM

# Preprocessor options
CPP = /usr/bin/cpp
FCPPFLAGS = $(ENVIRON)
FSRCDIR = ./srctmp

#directories
OD = ./obj
INCDIR = -I/usr/include -I$(BASE_NETCDF)/include 
LIBDIR = -L. -L$(BASE_NETCDF)/lib 

all: postproc 

QPSO: $(OD)/QPSO.o $(OD)/init_random_seed.o
	$(F90) -o $@ $(FFLAGS) $(INCDIR) $^ $(LIBDIR) -lnetcdf

ERGO_CLM: $(OD)/main.o $(OD)/ergo_global.o $(OD)/ERGO.o $(OD)/ergo_fileIO.o $(OD)/ergo_input_vars.o $(OD)/ergo_converge.o $(OD)/ergo_timer.o $(OD)/ergo_rndseed.o $(OD)/DE.o $(OD)/ceo.o $(OD)/mt19937ar.o $(OD)/model.o 
	$(F90) -o $@ $(FFLAGS) $(INCDIR) $^ $(LIBDIR) -lnetcdf


postproc: $(OD)/postproc.o
	$(F90) -o $@ $(FFLAGS) $(INCDIR) $^ $(LIBDIR) -lnetcdf

#object files
$(OD)/postproc.o: postproc.f90
	$(F90) -c $(FFLAGS) $(INCDIR) postproc.f90 -o $@
$(OD)/QPSO.o: QPSO_src/qpso.f90
	$(F90) -c $(FFLAGS) $(INCDIR) QPSO_src/qpso.f90 -o $@
$(OD)/init_random_seed.o: QPSO_src/init_random_seed.f90
	$(F90) -c $(FFLAGS) $(INCDIR) QPSO_src/init_random_seed.f90 -o $@

#ERGO optimization code
$(OD)/main.o: $(OD)/ergo_global.o $(OD)/ERGO.o $(OD)/ergo_fileIO.o $(OD)/model.o $(FSRCDIR)/main.f90
	$(F90) -c $(FFLAGS) $(INCDIR) $(FSRCDIR)/main.f90 -o $@
$(OD)/ergo_global.o: ERGO_src/ergo_global.f90
	$(F90) -c $(FFLAGS) $(INCDIR) ERGO_src/ergo_global.f90 -o $@
$(OD)/ERGO.o: $(OD)/ergo_global.o $(OD)/ergo_fileIO.o $(OD)/ergo_input_vars.o $(OD)/ergo_converge.o $(OD)/ergo_timer.o $(OD)/model.o $(OD)/DE.o $(OD)/ceo.o $(OD)/ergo_rndseed.o $(FSRCDIR)/ERGO.f90
	$(F90) -c $(FFLAGS) $(INCDIR) $(FSRCDIR)/ERGO.f90 -o $@
$(OD)/DE.o: $(OD)/ergo_global.o $(OD)/model.o $(OD)/ergo_converge.o $(OD)/ergo_timer.o $(OD)/mt19937ar.o $(OD)/ergo_rndseed.o $(FSRCDIR)/DE.f90	
	$(F90) -c $(FFLAGS) $(INCDIR) $(FSRCDIR)/DE.f90 -o $@
$(OD)/ceo.o: $(OD)/ergo_global.o $(OD)/model.o $(OD)/ergo_converge.o $(OD)/ergo_timer.o $(OD)/ergo_rndseed.o $(OD)/mt19937ar.o $(FSRCDIR)/ceo.f90
	$(F90) -c $(FFLAGS) $(INCDIR) $(FSRCDIR)/ceo.f90 -o $@
$(OD)/ergo_fileIO.o: $(OD)/ergo_global.o $(OD)/ergo_timer.o $(OD)/ergo_input_vars.o $(OD)/model.o $(FSRCDIR)/ergo_fileIO.f90
	$(F90) -c $(FFLAGS) $(INCDIR) $(FSRCDIR)/ergo_fileIO.f90 -o $@
$(OD)/ergo_converge.o: $(OD)/ergo_global.o $(OD)/model.o $(OD)/ergo_fileIO.o $(OD)/ergo_input_vars.o ERGO_src/ergo_converge.f90
	$(F90) -c $(FFLAGS) $(INCDIR) ERGO_src/ergo_converge.f90 -o $@
$(OD)/ergo_input_vars.o: $(OD)/ergo_global.o $(FSRCDIR)/ergo_input_vars.f90
	$(F90) -c $(FFLAGS) $(INCDIR) $(FSRCDIR)/ergo_input_vars.f90 -o $@
$(OD)/ergo_timer.o: $(OD)/ergo_global.o $(FSRCDIR)/ergo_timer.f90
	$(F90) -c $(FFLAGS) $(INCDIR) $(FSRCDIR)/ergo_timer.f90 -o $@
$(OD)/ergo_rndseed.o:  $(OD)/ergo_global.o $(OD)/mt19937ar.o $(FSRCDIR)/ergo_rndseed.f90
	$(F90) -c $(FFLAGS) $(INCDIR) $(FSRCDIR)/ergo_rndseed.f90 -o $@
$(OD)/mt19937ar.o:  ERGO_src/mt19937ar.f90
	$(F90) -c $(FFLAGS) $(INCDIR) ERGO_src/mt19937ar.f90 -o $@
$(OD)/model.o: $(OD)/ergo_rndseed.o $(OD)/ergo_global.o ERGO_src/model_CLM.f90
	$(F90) -c $(FFLAGS) $(INCDIR) ERGO_src/model_CLM.f90 -o $@
$(OD)/getposterior.o: ERGO_src/getposterior.f90
	$(F90) -c $(FFLAGS) $(INCDIR) ERGO_src/getposterior.f90 -o $@

# Preprocessing objects:
$(FSRCDIR)/main.f90 : ERGO_src/main.f90
	$(CPP) $(FCPPFLAGS) ERGO_src/main.f90 > $(FSRCDIR)/main.f90
$(FSRCDIR)/ERGO.f90 : ERGO_src/ERGO.f90
	$(CPP) $(FCPPFLAGS) ERGO_src/ERGO.f90 > $(FSRCDIR)/ERGO.f90
$(FSRCDIR)/ceo.f90 : ERGO_src/ceo.f90
	$(CPP) $(FCPPFLAGS) ERGO_src/ceo.f90 > $(FSRCDIR)/ceo.f90
$(FSRCDIR)/DE.f90 : ERGO_src/DE.f90
	$(CPP) $(FCPPFLAGS) ERGO_src/DE.f90 > $(FSRCDIR)/DE.f90
$(FSRCDIR)/ergo_fileIO.f90 : ERGO_src/ergo_fileIO.f90
	$(CPP) $(FCPPFLAGS) ERGO_src/ergo_fileIO.f90 > $(FSRCDIR)/ergo_fileIO.f90
$(FSRCDIR)/ergo_rndseed.f90 : ERGO_src/ergo_rndseed.f90
	$(CPP) $(FCPPFLAGS) ERGO_src/ergo_rndseed.f90 > $(FSRCDIR)/ergo_rndseed.f90
$(FSRCDIR)/ergo_input_vars.f90 : ERGO_src/ergo_input_vars.f90
	$(CPP) $(FCPPFLAGS) ERGO_src/ergo_input_vars.f90 > $(FSRCDIR)/ergo_input_vars.f90
$(FSRCDIR)/ergo_timer.f90 : ERGO_src/ergo_timer.f90
	$(CPP) $(FCPPFLAGS) ERGO_src/ergo_timer.f90 > $(FSRCDIR)/ergo_timer.f90


clean:
	rm -rf *.o *.mod obj/* postproc psens

depend:
	makedepend -- $(CFLAGS) -- $(SRCS) -- $(INCDIR)
