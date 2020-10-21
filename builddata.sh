#!/bin/bash
. ./findpdbabs.conf

if [ ! -d $tplsubdir ]; then
    mkdir $tplsubdir
fi;

(cd maketemplates; ./getpdbabseqs.pl)
(cd maketemplates; ./getgbtcrs.pl)
cat maketemplates/$tcrseqsfile $tplsubdir/*.faa > $tcrabsfile
\rm -f maketemplates/$tcrseqsfile

