#!/usr/bin/perl
use strict;
use lib '..';
use config;

my %config = config::ReadConfig("../findpdbabs.conf");

my $allseq="data/genbanktcrs.faa";
my $tmpseq="tmptcrseqs.faa";
my $repseq=$config{'tcrseqsfile'};

unlink $repseq;

# Grab and compile CD-HIT if not there
if( ! -d "cdhit" )
{
    print STDERR "Building cd-hit\n";
    `git clone git\@github.com:weizhongli/cdhit.git`;
    `(cd cdhit; make)`;
    print STDERR "done\n";
}

print STDERR "Running CD-HIT...";
`cdhit/cd-hit -c 0.6 -n 3 -i $allseq -o $tmpseq`;
print STDERR "done\n";

RelabelTCRs($repseq, $tmpseq);
unlink $tmpseq;
unlink "${tmpseq}.clstr";

sub RelabelTCRs
{
    my($repseq, $tmpseq) = @_;

    if(open(my $in, '<', $tmpseq))
    {
        if(open(my $out, '>', $repseq))
        {
            my $header   = '';
            my $sequence = '';
            while(<$in>)
            {
                chomp;
                if(/\>/)
                {
                    if($header ne '')
                    {
                        PrintFaa($out, $header, $sequence, 80);
                    }
                    s/\s.*//;                        # Remove from first space
                    $header = ">tcr" . substr($_,1); # Put tcr on the start
                    $sequence = '';
                }
                else
                {
                    $sequence .= $_;
                }
            }
            
            if($header ne '')
            {
                PrintFaa($out, $header, $sequence, 80);
            }
            close $out;
        }
        close $in;
    }
}

sub PrintFaa
{
    my($out, $header, $sequence, $minlen) = @_;

    if(length($sequence) >= $minlen)
    {
        print $out "$header\n";
        while(length($sequence))
        {
            my $this = substr($sequence, 0, 80);
            print $out "$this\n";
            $sequence = substr($sequence, 80);
        }
    }
}
