#!/bin/bash

 for i in `ls *`;do
    echo $i
    diff $i /glade/scratch/shpark/imsi/cesm/cesm2_0_beta04.4/components/cam/src/dynamics/mpas/external/$i
#    diff $i /glade/p/work/shpark/cesm/porting/201603/cesm1_4_beta07/components/cam/src/dynamics/mpas/external/$i
 done
