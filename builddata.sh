#!/bin/bash
. ./findpdbabs.conf

if [ ! -d $abtpldir ]; then
    mkdir -p $abtpldir
    if [ ! -d $abtpldir ]; then
        echo "No permision to create $abtpldir"
        exit 1
    fi
    echo "Created antibody template directory $abtpldir"
fi

mkdir -p $blastdir
mkdir -p $dbmdir

(cd maketemplates; ./getpdbabseqs.pl)
(cd maketemplates; ./getgbtcrs.pl)
cat maketemplates/$tcrseqsfile $tplsubdir/*.faa > $tcrabsfile
\rm -f maketemplates/$tcrseqsfile

