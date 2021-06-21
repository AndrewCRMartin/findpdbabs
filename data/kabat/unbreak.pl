#!/usr/bin/perl
use strict;

my $file = shift @ARGV;

my $labelSize = FindWidestLabel($file);

if(open(my $in, '<', $file))
{
    my $firstLine = 1;
    while(<$in>)
    {
        chomp;
        if(/^>(.*?)\|/)
        {
            my $id = $1;
            if(!$firstLine)
            {
                print "\n";
            }
            else
            {
                $firstLine = 0;
            }
            print "$id";
            print ' ' x (2 + $labelSize - length($id));
        }
        else
        {
            print;
        }
    }
    print "\n";
}

sub FindWidestLabel
{
    my($file) = @_;
    my $maxLabel = 0;

    if(open(my $in, '<', $file))
    {
        while(<$in>)
        {
            chomp;
            s/\s+$//;
            if(/^>(.*)/)
            {
                $maxLabel = length($1) if(length($1) > $maxLabel);
            }
        }
        close $in;
    }
    return($maxLabel);
}
