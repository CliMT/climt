#!/bin/bash
# To check if the actual COSP code on github differs from this version run
# check_cosp_github.sh and then see if there are any differences by having 
# a look at cosp_diff_to_github.out

# Get latest github version and unpack it
rm -rf COSPv2.0-master
rm -f master.zip
wget -nc https://github.com/CFMIP/COSPv2.0/archive/refs/heads/master.zip
if [ -f master.zip ] ; then
  unzip master.zip
  cp COSPv2.0-master/driver/data/inputs/UKMO/cosp_input_um_2d.nc .
else
  echo 'COSP github not accessible'
  exit 0
fi

# Now check the diff to the current fcm version
ierr=0
rm -f cosp_diff_to_github.out
diff -r -x cosp2_test.f90 COSPv2.0-master/driver/src ${RAD_DIR}/src/cosp_github/driver/src >> cosp_diff_to_github.out || ierr=1
diff -r -x README COSPv2.0-master/model-interface ${RAD_DIR}/src/cosp_github/model-interface >> cosp_diff_to_github.out || ierr=1
diff -r -x cosp_rttov_interface.F90 -x cosp_rttov.F90 -x cosp_rttov11.F90 COSPv2.0-master/src ${RAD_DIR}/src/cosp_github/src >> cosp_diff_to_github.out || ierr=1
diff -r COSPv2.0-master/subsample_and_optics_example ${RAD_DIR}/src/cosp_github/subsample_and_optics_example >> cosp_diff_to_github.out || ierr=1

if [ $ierr -gt 0 ] ; then
  echo 'COSP code differs from github version'
  exit 1
else
  rm -f cosp_diff_to_github.out
  rm -rf COSPv2.0-master
  rm -f master.zip
  echo 'COSP code matches github version'
  exit 0
fi
