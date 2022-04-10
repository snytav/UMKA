FC=mpif90
CFLAGS=-I.
DEPS = mpif.h


%.o: %.f $(DEPS)
	$(FC) -c -o $@ $< $(CFLAGS)

all:  beam.o plasma.o para.o vvod.o output.o energy.o fourier.o out3d.o \
     check.o control.o micro.o
	$(FC) -o th *.o 
	
clean:
	rm *.o th
