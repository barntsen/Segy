#!/bin/sh
# Run tests of segy/su read/write routines 
cp ../Src/segy.py .
cp ../Src/ibm_float.py .

#echo "=== segyread ======="
#python3 segyread.py

echo "=== suread ======="
python3 suread.py

