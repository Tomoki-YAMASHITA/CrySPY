#note the option for the fort compiler
#F90=ifort
F90=gfortran
F90_OPTS= -static-libgfortran  -cpp -DUSE_M_CMD_ARG

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------
SRC = precmod.f90 stringmod.f90  constant.f90 \
      findinv.f90 convolute.f90   cal_fingerprint.f90   

OBJ = $(SRC:.f90=.o)

all:  $(OBJ)  cal_fingerprint 


cal_fingerprint:	$(OBJ)
	$(F90) $(F90_OPTS) -o cal_fingerprint cal_fingerprint.o  precmod.o stringmod.o constant.o convolute.o  findinv.o

cleanall:
	@rm -f *.o *.mod gen_kmesh  rndgen 
clean:
	@rm -f *.o *.mod 

#end
