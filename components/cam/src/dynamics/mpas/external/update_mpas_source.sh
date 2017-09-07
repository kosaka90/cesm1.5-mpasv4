#!/usr/bin/env sh

################################################################################
#
# This script may be used to obtain updated MPAS non-hydrostatic atmosphere
#    code from a compiled working copy of the MPAS repository code.
# The environment variable MPAS_PATH should be set to the path of the top-level 
#    MPAS directory.
#
# Initial version: 6 June 2010
#
################################################################################


MPAS_PATH=${MPAS_PATH:-/glade/p/work/duda/cam_mpas_nh}

if [ ! -d $MPAS_PATH ]; then
   echo "Specified MPAS_PATH not found: $MPAS_PATH"
   exit 1
fi

if [ ! -e $MPAS_PATH/src/driver/mpas_cam_interface.f90 ]; then
   echo "MPAS must be compiled before running this script"
   exit 1
fi

rm -f *.c *.h *.f90

cp ${MPAS_PATH}/src/driver/mpas_cam_interface.f90 .
cp ${MPAS_PATH}/src/framework/*.c .
cp ${MPAS_PATH}/src/external/ezxml/*.c .
cp ${MPAS_PATH}/src/external/ezxml/*.h .
cp ${MPAS_PATH}/src/framework/*.f90 .
cp ${MPAS_PATH}/src/operators/*.f90 .
cp ${MPAS_PATH}/src/core_atmosphere/*.f90 .
cp ${MPAS_PATH}/src/core_atmosphere/dynamics/*.f90 .

#mv mpas_atm_advection.f90 atm_advection.f90
#mv mpas_atm_mpas_core.f90 mpas_core.f90
#mv mpas_atm_time_integration.f90 atm_time_integration.f90

rm -f mpas_tracer_advection_*

