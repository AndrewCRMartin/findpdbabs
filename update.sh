#!/bin/bash

thedate=`date +%Y-%m-%d`
echo "Extracting sequences for new PDB files..."
make -f Makefile.pdbseq

echo "Identifying antibodies"
./findpdbabs.pl newabs_${thedate}.out
cat newabs_${thedate}.out >>abs.out 
