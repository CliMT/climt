#!/bin/bash
# Script to run both LW and SW spectral file generation for TRAPPIST-1e

set -e

echo "Starting TRAPPIST-1e spectral file generation pipeline..."

# Run LW script
echo "--------------------------------------------------------"
echo "Running Longwave (LW) generation..."
./mk_sp_lw_17_trappist1e

# Run SW script
echo "--------------------------------------------------------"
echo "Running Shortwave (SW) generation..."
./mk_sp_sw_42_trappist1e

echo "--------------------------------------------------------"
echo "All spectral files generated successfully."
