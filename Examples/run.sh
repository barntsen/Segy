#!/bin/sh
# Run tests of segy/su read/write routines 
cp ../Src/segy.py .
cp ../Src/ibm_float.py .

echo "=== suwrite ======="
python3 suwrite.py 
echo "=== suread ======="
python3 suread.py 

echo "=== segywrite ======="
python3 segywrite.py
echo "=== segyread ======="
python3 segyread.py

