#!/bin/bash

if [ ! -f ../findpdbabs.conf ]; then
    echo "You must create a configuration file: findpdbabs.conf"
    exit 1
fi

if [ $0 != "./builddata.sh" ]; then
    echo "Script must be run in the dataprep directory as ./builddata.sh"
    exit 1
fi

. ../findpdbabs.conf

function makedir
{
    dir=$1
    msg=$2
    if [ ! -d $dir ]; then
        mkdir -p $dir
        if [ ! -d $dir ]; then
            echo "No permision to create $dir"
            exit 1
        fi
        echo "Created $msg ($dir)"
    else
        echo "$msg already exists ($dir)"
    fi
}

makedir $abtpldir "Antibody template directory"
makedir $blastdir "Blast database directory"
makedir $dbmdir   "DBM directory"

(cd maketemplates; ./getpdbabseqs.pl)
(cd maketemplates; ./getgbtcrs.pl)
cat $tcrseqsfile $abtpldir/*.faa > $tcrabsreffile
#\rm -f $tcrseqsfileq

