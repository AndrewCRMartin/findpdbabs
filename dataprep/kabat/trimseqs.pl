#!/usr/bin/perl

use strict;
my $skipStart = shift(@ARGV);
my $skipEnd   = shift(@ARGV);

while(<>)
{
    chomp;
    
    my($id, $seq) = split;

    $seq = substr($seq, $skipStart, ($skipEnd-$skipStart));

    printf "%-25s", $id;
    print "$seq\n";
}
