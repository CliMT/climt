
# IF YOU ARE REALLY SURE that you want to use the github version rather
# than the current fcm version, you can use copy_cosp_github.sh
# to overwrite relevant files with their github versions.

# Get latest github version and unpack it
rm -rf COSPv2.0-master
rm -f master.zip
wget -nc https://github.com/CFMIP/COSPv2.0/archive/refs/heads/master.zip
unzip master.zip

# Now overwrite with the github version
cp -R COSPv2.0-master/driver/src ${RAD_DIR}/src/cosp_github/driver
cp -R COSPv2.0-master/model-interface ${RAD_DIR}/src/cosp_github
cp -R COSPv2.0-master/src ${RAD_DIR}/src/cosp_github
cp -R COSPv2.0-master/subsample_and_optics_example ${RAD_DIR}/src/cosp_github
