#!/usr/bin/perl

use strict;
use fasta;

my $inFile          = shift @ARGV;
my $db              = shift @ARGV;
my %datafiles       = ();
$datafiles{'alpha'} = shift @ARGV;
$datafiles{'beta'}  = shift @ARGV;
$datafiles{'gamma'} = shift @ARGV;
$datafiles{'delta'} = shift @ARGV;

my $fasta="/usr/local/apps/fasta33/fasta33_t";

if(open(my $in, '<', $inFile))
{
    
    my($id, $info, $sequence);

    while((($id, $info, $sequence) = fasta::ReadFasta($in)) && ($id ne ''))
    {
        my $tfile = "/tmp/rfa_" . $$ . ".faa";
        if(open(my $out, '>', $tfile))
        {
            fasta::PrintFasta($out, $id, $sequence);
            close $out;
            my $result = `$fasta -q -p -b 10 -d 0 -E 15 $tfile $db 2`;
            my $dataset = FindDataset($result);

            if(defined($datafiles{$dataset}))
            {
                if(open(my $fp, '>>', $datafiles{$dataset}))
                {
                    fasta::PrintFasta($fp, "$id|$dataset", $sequence);
                    close $fp;
                }
            }
        }
    }
    close $in;
}

sub FindDataset
{
    my($result) = @_;

    my @lines = split(/\n/, $result);
    my $inData = 0;
    foreach my $line (@lines)
    {
        if($line =~ /The best scores are/)
        {
            $inData = 1;
        }
        elsif($inData)
        {
            my(@fields) = split(/[\|\s+]/, $line);
            return($fields[1]);
        }
    }
    return('unknown');
}
