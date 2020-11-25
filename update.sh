#!/bin/bash

if [ ! -f Makefile.pdbseq ]; then
    ./configure.pl
fi

abfile=abs.out

thedate=`date +%Y-%m-%d`
echo "*** Extracting sequences for new PDB files..."
make -f Makefile.pdbseq

echo ""
echo "*** Identifying antibodies..."
newabfile=newabs_${thedate}.out
./findpdbabs.pl $newabfile
cat $newabfile >>$abfile

echo ""
echo "********************************************************************"
echo "New antibodies are listed in $newabfile"
echo "All antibodies are listed in $abfile"
echo "********************************************************************"

