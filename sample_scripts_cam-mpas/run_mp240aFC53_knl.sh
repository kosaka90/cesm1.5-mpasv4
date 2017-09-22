#!/bin/bash
#run this script in $CESM_ROOT/cime/scripts. After creating the case, 
#I usually create a copy under the case directory, and run from there.

casename="test01_240aFC53_knl"

CESM_ROOT="/global/u1/k/ksa/MPAS/model/cesm2_0_beta05"

SIM_ROOT="/global/cscratch1/sd/ksa/simulation/MPAS"

CASE_ROOT="/global/project/projectdirs/m1867/ksa/MPAS/cases"

CASE_RES="mp240a_g16"

CASE_MACH="cori-knl"

MPAS_INPUT_DIR="/global/project/projectdirs/m1867/MPASinput/mp240a"

CAM_VER="cam5"

AERO_OPT="none"

echo "CESM_ROOT=" $CESM_ROOT
echo "SIM_ROOT=" $SIM_ROOT
echo "CASE_ROOT=" $CASE_ROOT
echo "CASE_RES=" $CASE_RES
echo "CASE_MACH=" $CASE_MACH
echo "MPAS_INPUT_DIR=" $MPAS_INPUT_DIR   #this is only within the script process under cshell terminal?
echo "CAM version = " $CAM_VER
echo "Aerosol scheme = " $AERO_PROG

#additionally set env.var for land model input 
MPAS_LNDINPUT_DIR=${MPAS_INPUT_DIR}"/clm4_0"
echo "MPAS_LNDINPUT_DIR=" $MPAS_LNDINPUT_DIR
MPAS_RUN_DIR="$SIM_ROOT/$casename/run"
echo "run directory is " $MPAS_RUN_DIR


#create a case:
if [ -d $CASE_ROOT/$casename ]
then
    echo "case directory exists"
else
    cd $CESM_ROOT/cime/scripts
    pwd
    echo "creating a new case"
    ./create_newcase -case $CASE_ROOT/$casename -compset FC5 -res $CASE_RES -mach $CASE_MACH
fi

cd $CASE_ROOT/$casename
#./create_newcase -case /global/project/projectdirs/m1867/ksa/MPAS/cases/test01_240aFC53_knl -compset FC5 -res mp240a_g16 -mach cori-knl

casedir=`pwd`
#use backtick ` to capture command output to a variable


ntasks=68   # number of MPI tasks. Note MPAS does not support openMP threads
ppnode=68   # number of logical cores (MPI ranks) to use per node (control hyperthreads) 

#env_mach_pes.xml. 
#Need to run case.setup (and --clean) everytime (except first) we edit this file
edit_pes=true
if [ "$edit_pes" = true ]; then
    echo "editing env_mach_pes.xml"
    $casedir/xmlchange NTASKS_ATM=$ntasks
    $casedir/xmlchange NTASKS_OCN=$ntasks
    $casedir/xmlchange NTASKS_LND=$ntasks
    $casedir/xmlchange NTASKS_ROF=$ntasks
    $casedir/xmlchange NTASKS_ICE=$ntasks
    $casedir/xmlchange NTASKS_GLC=$ntasks
    $casedir/xmlchange NTASKS_WAV=$ntasks
    $casedir/xmlchange NTASKS_CPL=$ntasks
    $casedir/xmlchange PES_PER_NODE=$ppnode        
fi
$casedir/xmlquery NTASKS
$casedir/xmlquery NTHRDS
$casedir/xmlquery PES_PER_NODE


#env_build.xml. 
#must run the case.build command after changing this file.
#safe to run clean_build before case.build
edit_build=true
if [ "$edit_build" = true ]; then
    echo "editing env_build.xml"
    $casedir/xmlchange CIME_OUTPUT_ROOT=$SIM_ROOT
    $casedir/xmlchange DEBUG=FALSE
    #$casedir/xmlchange CAM_DYCORE=fv  MPAS is set as the default dycore for this CESM package
    $casedir/xmlchange -id CLM_CONFIG_OPTS -val ' -phys clm4_0'
    $casedir/xmlchange -id CAM_CONFIG_OPTS -val " -phys $CAM_VER  -chem $AERO_OPT"

fi
$casedir/xmlquery DEBUG
$casedir/xmlquery BUILD_THREADED
$casedir/xmlquery SMP_BUILD
$casedir/xmlquery CAM_DYCORE
$casedir/xmlquery CAM_CONFIG_OPTS

#env_run.xml
#Sets run-time settings such as length of run, frequency of restarts, output of coupler diagnostics,
#and short-term and long-term archiving. This file can be edited at any time before a job starts.
edit_run=true
if [ "$edit_run" = true ]; then
    echo "editing env_run.xml"
    $casedir/xmlchange RUN_STARTDATE=2000-01-01
    $casedir/xmlchange STOP_OPTION=ndays
    $casedir/xmlchange STOP_N=5
    $casedir/xmlchange DOUT_S=FALSE
    $casedir/xmlchange INFO_DBUG=1
    #$casedir/xmlchange ATM_NCPL=144   #set physics time step by this variable
    #$casedir/xmlchange REST_OPTION='$STOP_OPTION'
    $casedir/xmlchange REST_OPTION=never
    #$casedir/xmlchange REST_N='$STOP_N'
    $casedir/xmlchange HIST_OPTION=never
fi
#check
$casedir/xmlquery RUN_STARTDATE
$casedir/xmlquery STOP_OPTION
$casedir/xmlquery STOP_N
$casedir/xmlquery DOUT_S



#env_batch.xml. 
#Sets batch system specific settings such as wallclock time and queue name.
#Can edit anytime
#but editing this way makes the queue and wall clock time same for all types of jobs
#run, short-term archive, and long-term archive? OK for testing, better to edit manually....
edit_batch=true
if [ "$edit_batch" = true ]; then
    echo "editing env_batch.xml"
    $casedir/xmlchange PROJECT=m1867
    $casedir/xmlchange JOB_WALLCLOCK_TIME=00:30:00
    $casedir/xmlchange JOB_QUEUE=debug
    #$casedir/xmlchange JOB_WALLCLOCK_TIME=02:00:00
    #$casedir/xmlchange JOB_QUEUE=regular

fi      
#check
$casedir/xmlquery JOB_WALLCLOCK_TIME
$casedir/xmlquery PROJECT
$casedir/xmlquery JOB_QUEUE

#env_mach_specific.xml
#Sets a number of machine-specific environment variables for building and/or running.
#The variables in this file do not have "id" and "value" attribute, so can't be edited 
#by xmlchange
#always check cpu_bind ("cores" when ntasks == ppnode otherwise "threads"), depending on pes setting


### namelist options
# CAM namelist. For MPAS, we need phys_loadbalance=0 
# MPAS initial condition does not have aerosols, so need state_debug_checks=.false.
# for prescribed aerosol (with MG2 only?), need  use_hetfrz_classnuc = .false.
# and other MPAS specific input files need to be specified
cat >| user_nl_cam <<HERE
phys_loadbalance=0
state_debug_checks=.false.
use_hetfrz_classnuc = .false.
empty_htapes=.true.
fincl1 = 'U:I','OMEGA:I','PRESSURE:I','Z3:I'
nhtfrq = -1
mfilt = 24
ncdata	='$MPAS_INPUT_DIR/x1.10242.CAM.init.161018.nc'
bnd_topo='$MPAS_INPUT_DIR/mp240a_topo_150210.nc'
drydep_srf_file='$MPAS_INPUT_DIR/atmsrf_mp240a_150519.nc'
HERE

# CLM namelist. CLM runs on the MPAS grid, so it requires regridded input data
cat >| user_nl_clm <<HERE
hist_empty_htapes=.true.
finidat='$MPAS_LNDINPUT_DIR/clmi.BCN.2000-01-01_mp240a_gx1v6_simyr2000_c150311.nc'
flanduse_timeseries='$MPAS_LNDINPUT_DIR/landuse.timeseries_mp240a_hist_simyr1850-2005_c150311.nc'
HERE

#CPL namelist
cat >| user_nl_cpl <<HERE
aoflux_grid='atm'
HERE

#### run case.setup ####
run_setup=true
clean_setup=false    #for changing PE layout
if [ "$run_setup" = true ]; then
    if [ "$clean_setup" = true ]; then
        ./case.setup --clean
    fi

    ./case.setup
fi


###copy modified source codes into the local SourceMod directory for each case
echo "copying modified source codes into $case/SourceMods/"
#Take care of the number of advected scalars
#for prescribed aerosols with MG1 (CAM5.3)
if [ $AERO_OPT = "none" ];then 
    if [ $CAM_VER = "cam5" ];then 
       echo "copying mpas_atm_core_interface.f90 for prescribed aerosol with MG1 (CAM5.3)"
       /usr/bin/cp -f /global/project/projectdirs/m1867/MPASinput/SourceMods/PM_MG1_MPASv4/src.cam/mpas_atm_core_interface.f90 $casedir/SourceMods/src.cam/

    else
       echo "copying mpas_atm_core_interface.f90 for prescribed aerosol with MG2 (CAM5.4 or later)"

       /usr/bin/cp -f /global/project/projectdirs/m1867/MPASinput/SourceMods/PM_MG2_MPASv4/src.cam/mpas_atm_core_interface.f90 $casedir/SourceMods/src.cam/
    fi
fi


#### run case.build ####
run_build=true
clean_build=false
if [ "$run_build" = true ]; then
    if [ "$clean_build" = true ]; then
        ./case.build --clean    #for modifying source codes 
        #./case.build --clean-all  #for changing PE layout
    fi

    ./case.build
fi

#### copy MPAS-specific namelist into the run directory
cp_mpasin=false
if [ "$cp_mpasin" = true ]; then
    cd ${MPAS_RUN_DIR}
    pwd
    ln -s ${MPAS_INPUT_DIR}/stream_list.atmosphere.diagnostics
    ln -s ${MPAS_INPUT_DIR}/stream_list.atmosphere.output
    ln -s ${MPAS_INPUT_DIR}/stream_list.atmosphere.surface
    cp ${MPAS_INPUT_DIR}/streams.atmosphere .
    cp ${MPAS_INPUT_DIR}/namelist.atmosphere .

fi

echo "make sure to check the MPAS input namelist namelist.atmosphere in ${MPAS_RUN_DIR}"
echo "use case.submit to submit a job"

