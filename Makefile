# Sample makefile

FC=gfortran.exe

all: 
	$(FC) $(cdebug) $(cflags) $(cvars) beam.f check.f control.f energy.f fourier.f out3d.f output.f para.f plasma.f vvod.f

clean:
	del a.exe

#simple.exe: simple.obj
#  $(link) $(ldebug) $(conflags) -out:simple.exe simple.obj $(conlibs) lsapi32.lib

#challeng.exe: challeng.obj md4c.obj
#  $(link) $(ldebug) $(conflags) -out:challeng.exe $** $(conlibs) lsapi32.lib