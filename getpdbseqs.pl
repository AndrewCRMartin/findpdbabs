#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    findpdbabs
#   File:       getpdbseqs.pl
#   
#   Version:    V1.0
#   Date:       25.11.20
#   Function:   Update sequence files for new PDB files
#   
#   Copyright:  (c) Prof. Andrew C. R. Martin, UCL, 2020
#   Author:     Prof. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#
#*************************************************************************
# Add the path of the executable to the library path
use FindBin;
use lib $FindBin::Bin;
use strict;
use config;

my $configFile = 'findpdbabs.conf';
my $exeDir     = $FindBin::Bin;
#*************************************************************************
my %config = config::ReadConfig($configFile);

if(opendir(my $dh, $config{'pdbdir'}))
{
    my @pdbfiles = grep(!/^\./, readdir($dh));
    closedir($dh);

    my $nConverted = 0;
    foreach my $pdbfile (@pdbfiles)
    {
        my $seqfile = $pdbfile;
        $seqfile =~ s/\.ent/.faa/;
        $seqfile = $config{'faadir'} . "/$seqfile";
        if(! -f $seqfile)
        {
            my $fullpdbfile = $config{'pdbdir'} . "/$pdbfile";
            print "Extracting sequence from $fullpdbfile to $seqfile\n";
            `pdb2pir -s -c -f $fullpdbfile > $seqfile`;
            $nConverted++;
        }
    }

    print "\n\n$nConverted sequence files were created from PDB files\n";
}
else
{
    print STDERR "Can't open directory $config{'pdbdir'}\n";
}
