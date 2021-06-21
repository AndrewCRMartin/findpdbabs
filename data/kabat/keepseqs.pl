#!/usr/bin/perl

use strict;
my $type      = shift(@ARGV);
my $skipStart = shift(@ARGV);
my $skipEnd   = shift(@ARGV);

while(<>)
{
    chomp;
    
    my($id, $seq) = split;
    if(! ($seq=~/X/))
    {
        my $startDash = 0;
        my $endDash   = 0;
        if($seq =~ /^(-+)/)
        {
            $startDash = length($1);
        }
        if($seq =~ /(-+)$/)
        {
            $endDash = length($1);
        }
        if(($startDash <= $skipStart) &&
           ($endDash   <= $skipEnd))
        {
            print ">$id|$type\n";
            $seq =~ s/\-//g;
            print "$seq\n";
        }
    }
}
