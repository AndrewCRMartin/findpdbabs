#!/bin/bash

thedate=`date +%Y-%m-%d`
make -f Makefile.pdbseq
./findpdbabs.pl >newabs_${thedate}.out 2>findpdbabs.log
cat newabs_${thedate}.out >>abs.out 
