RunDir=test          # WORK DIRECTORY. 

# GENERAL SETTINGS
proc=1                   # NUMBER OF PROCESSORS FOR PARTICLES
procf=1                  # TOTAL NUMBER OF PROCESSORS
WTIME=144                # WORKTIME IN HOURS
NL=100                   # LONGITUDINAL NUMBER OF NODES
NY=4                     # TRANSVERSE NUMBER OF NODES: Y
NZ=4                     # TRANSVERSE NUMBER OF NODES: Z
nt1=1                   # INNER ITERATIONS
nt2=1     #dnee         # OUTER ITERATIONS
nt3=5                    # nt2/nt3 beam                
ntd=10000000             # MOMENT TO START MORE DETAILED OUTPUT
ntd1=40                  # DETAILED nt1 value 
tau=0.001d0                # TIMESTEP
var=0                    # TYPE number for parameter variation
magf=1                   # MAGNETIC FIELD ON/OFF
beamf=1                  # BEAM ON/OFF 
elmod=0e0                # ELECTRIC FIELD MODULATION ON/OFF
smooth=0                 # BEAM CURRENT SMOOTH MODE (1 - 5-point smooth, 2 - averaging)

# PHYSICAL PARAMETERS
xm=1.1424                # X DOMAIN SIZE
ym=0.05                  # Y DOMAIN SIZE
zm=0.05                  # Z DOMAIN SIZE
xmp=$xm                  # X PLASMA SIZE
ymp=0.025                 # Y PLASMA SIZE
zmp=0.05                 # Z PLASMA SIZE
xmb=0.2                  # X BEAM SIZE
ymb=0.025                 # Y BEAM SIZE
zmb=0.05                 # Z BEAM SIZE
rbd=2e-4                 #  BEAM TO ELECTRON DENSITY RATIO 
rimp=0.2                 # BEAM IMPULSE negative - NO BEAM
xr1=0.1d0                # ELECTRON TEMPERATURE 
alp=1.0d-2               # PLASMA DENSITY
lp=10                    # PARTICLES PER CELL (EACH SORT)
mod=0d0                  # FREQUENCY OF INITIAL VELOCITY MODULATION
melb=1.0e0               # MASS OF BEAM ELECTRON
mel=1.0e0                # RELATIVE MASSES OF MODEL ELECTRONS
mi=1836e0                # AND IONS 
ni0=1.0e0                # INITIAL ION DENSITY                 1 - 1e14 cm-3
hx0=0.0                  # INITIAL X MAGNETIC FIELD
Tex0=1e-3                # INITIAL X ELECTRON TEMPERATURE        1 - 560 eV
Tey0=$Tex0               # INITIAL Y ELECTRON TEMPERATURE        1 - 560 eV
Tez0=$Tex0               # INITIAL Z ELECTRON TEMPERATURE        1 - 560 eV
Tbeam=0.14               # BEAM VELOCITY DELTA/beam velocity (tablichnoe znachenie)   =delta

#OUTPUT FLAGS
fte=0                         # ELECTRON TEMPERATURE FOURIER
fti=0                         # ION TEMPERATURE FOURIER
fne=0                         # ELECTRON DENSITY FOURIER
fnb=0                         # BEAM ELECTRON DENSITY
fee=0                         # ELECTRIC FIELD FOURIER
fni=0                         # ION DENSITY FOURIER
fhh=0                         # MAGNETIC FIELD FOURIER
map=0                         # ALL THE ABOVE LISTED AS SMALL SIZE MAPS
rdm=0                         # DEBYE RADIUS MAP
rlm=0                         # LARMOR RADIUS MAP
radf=0                        # RADIAL AVERAGE DISTRIBUTIONS
outterm=0                     # PRINT EVERY TIMESTEP: TEMPERATURE AND ENERGIES
outf=0                        # FOURIER OUTPUT ON/OFF
ons=0                         # GENERAL FOURIER AND TEMPERATURE NUMBERS REPORTED ONLY ONCE
disf=0                        # 1D DISTRIBUTION FUNCTION
o3d=0                         # 3D heat flow
ChainOut=0                    # CURRENTS AND ELECTRIC FIELD
curout=0                      # X CURRENT
EnOut=0                       # VARIOUS KINDS OF ENERGY
CurBal=0                      # BALANCE OF CURRENTS (BEAM, ELECTRON AND ION)
Dens1D=0                      # DENSITIES ALONG X
ctrl_attr=0                   # WRITE CONTROL ATTRIBUTES FILE
write_all=0                   # WRITE ALL VALUES IN A BINARY FILE
beam_write=0                  # WRITE BEAM VALUES
jx_four=0                     # X Current Fourier Transform
wpbeam=0                      # WRITE ALL BEAM PARTICLES
wpelec=0                      # WRITE ALL ELECTRON PARTICLES
wpion=0                       # WRITE ALL ION PARTICLES
baddst=0                      # EXTRA BEAM PARTICLES STATISTICS
 

# UNKNOWN PARAMETERS
hx0=0.0                       
pak=1 
sst=1.0                        
u0=0.0           
dT=0.1 
kk=0.001                      

# TECHNICAL PARAMETERS
trpart=5                # PORTION OF PARTICLES TRANSMITTED AT ONE TIMESTEP
deb=0                   # DEBUG LEVEL
jdeb=0                  # PROBE PARTICLE
omp=100                   # OUTPUT MAP SIZE
growtot=10              # TOTAL NUMBER OF MOSTLY GROWING HARMONICS BEING TRACED   


# END OF USER-EDITED PART
#*****************************************************************************
#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


echo "      integer nproc,fproc
            parameter(nproc="$proc",fproc="$procf")">&procs.par
echo "      integer imp
            parameter(imp="$[$NL+2]")">&laser1.par
echo "      include 'procs.par'
            integer lmp
            parameter(lmp="$NY"/nproc+2)">&laser2.par
echo "      integer kmp
            parameter(kmp="$[$NY+2]")">&laser3.par
echo "      integer jmp
            parameter(jmp="$[2*$NL*$NY*$NZ*$lp]"/fproc)">&laser4.par
echo "      real*8 mel,mi
            parameter(mel="$mel",mi="$mi") ">&mass.par
	    
echo $var>&var.dat	    

echo "      include 'laser1.par'
            include 'laser2.par'
            include 'laser3.par'
            include 'laser4.par'
   	    integer deb,lmp2,bjmp,jdeb,outf,lmpf,ons,omp,outterm
	    integer fte,fti,fee,fhh,map,rdm,rlm,fne,fni,m7nt2,attr
	    integer disf,beamf,magf,grtot,out3d,smooth,ctrl_attr
	    integer radf,totsteps,ntd,ntd1,lnt1,fnb,outfsrc,write_all
	    real*8 beta0,tex0,tey0,tez0,rbd,rimp,lx0,ly0,lz0,Tb,mfrq
	    real*8 melb,elmod
            real*8 xmp,ymp,zmp 
            real*8 xmb,ymb,zmb 
	    integer jx_four
	    integer chnout,enout,curbal,out1dd,curout,beam_write
	    integer wpbeam,wpelec,wpion,baddst
	    parameter(wpbeam="$wpbeam",wpelec="$wpelec",wpion="$wpion")
            parameter(xmp="$xmp",ymp="$ymp",zmp="$zmp",baddst="$baddst")
            parameter(xmb="$xmb",ymb="$ymb",zmb="$zmb")
	    parameter(out1dd="$Dens1D",curout="$curout")
	    parameter(enout="$EnOut",curbal="$CurBal",write_all="$write_all")
	    parameter(outfsrc=1,chnout="$ChainOut",ctrl_attr="$ctrl_attr")
	    parameter(melb="$melb",elmod="$elmod",smooth="$smooth")
	    parameter(lx0="$xm",Tb="$Tbeam",mfrq="$mod")
	    parameter(ly0="$ym",lz0="$zm",ntd="$ntd",ntd1="$ntd1")
	    parameter(beta0=1d0,rbd="$rbd*$beamf",rimp="$rimp",grtot="$growtot")
	    parameter(beamf="$beamf",magf="$magf",lnt1="$nt1")
	    parameter(tex0="$Tex0",tey0="$Tey0",tez0="$Tez0")
            parameter(bjmp=jmp/"$trpart",totSteps="$nt1*$nt2",fnb="$fnb")
	    parameter(fne="$fne",fni="$fni",radf="$radf",disf="$disf")
	    parameter(deb="$deb",outf="$outf",omp="$omp",ons="$ons",jdeb="$jdeb")
	    parameter(fte="$fte",fti="$fti",fee="$fee",fhh="$fhh")
	    parameter(map="$map",rdm="$rdm",rlm="$rlm",out3d="$o3d")
	    parameter(outterm="$outterm",m7nt2="$nt2",beam_write="$beam_write")
            parameter(lmp2=nproc*(lmp-2)+1,lmpf=nproc*(lmp-2)+2)
            parameter(jx_four="$jx_four")
	    integer me,mev,meh,vertComm,horComm
	    common/globme/me,mev,meh,vertComm,horComm ">&part.pf



# PHYSICAL PARAMETERS FILE
echo $nt1"                   nt1
"$nt2"                   nt2   
"$nt3"                   nt3   
0                        ml
"$tau"                   tau
"$lp"                    lp
"$ni0"                   ni0 
"$Te0"                   Te0
"$hx0"                   hx0
"$xm"                    lx0
"$ym"                    ly0 
"$zm"                    ly0 
"$pak"                   pak
"$sst"                    sst
"$u0"                    u0
"$dT"                    dT
"$kk"                    kk ">&lstart18.dat

    WorkDir=$RunDir'_mi_'$mi'_mel_'$mel'_melb_'$melb'_mfrq_'$mod'_Tex0_'$Tex0'_rbd_'$rbd'_rimp_'$rimp'_NL_'$NL'_NY_'$NY'_NZ_'$NZ'_TB_'$Tbeam'_xm_'$xm'_ym_'$ym'_zm_'$zm'_xmp_'$xmp'_ymp_'$ymp'_zmp_'$zmp'_xmb_'$xmb'_ymb_'$ymb'_zmb_'$zmb'_ni0_'$ni0'_tau_'$tau'_lp_'$lp'_procf_'$procf'_proc_'$proc'_bm_'$beamf
#     WorkDir=$RunDir'_mfrq_'$mod'_melb_'$melb'_NL_'$NL'_tau_'$tau'_lp_'$lp'_procf_'$procf'_proc_'$proc'_bm_'$beamf

    
     rm -r $WorkDir
     mkdir $WorkDir
    
    # COPYING FILES
    cp lstart18.dat ./$WorkDir
    cp var.dat ./$WorkDir
    cp *.par ./$WorkDir
    cp *f ./$WorkDir

    icc -c -g micro.c
    mpiifort -g -c *.f 2>&1 >& err.out
    mpiifort -g -o th *.o
    cat err.out
#    rm *.par *.dat
    
    cp th ./$WorkDir
    ./wnuc $procf $WTIME
    cp nuc ./$WorkDir
    cp zuf4 ./$WorkDir
    cd ./$WorkDir
     mpirun -np $procf ./th
#    qsub nuc
