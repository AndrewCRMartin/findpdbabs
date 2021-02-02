#!/usr/bin/perl
use strict;
use lib '../..';
use config;

my %config = config::ReadConfig("../../findpdbabs.conf");

my $allseq=shift (@ARGV);
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

RelabelTCRs($tmpseq, $allseq);

print STDERR "Running CD-HIT...";
`cdhit/cd-hit -T 0 -c 0.6 -n 3 -i $tmpseq -o $repseq`;
print STDERR "done\n";

unlink $tmpseq;
unlink "${repseq}.clstr";

sub RelabelTCRs
{
    my($repseq, $tmpseq) = @_;

    if(open(my $in, '<', $tmpseq))
    {
        if(open(my $out, '>', $repseq))
        {
            my $header     = '';
            my $sequence   = '';
            while(<$in>)
            {
                chomp;
                if(/\>/)
                {
                    if($header ne '')
                    {
                        PrintFaa($out, $header, $sequence, 80);
                    }

                    $header   = $_;
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

sub FixHeader
{
    my($header) = @_;
    $header = substr($header,1);
    $header =~ s/\s.*//;    # Remove from first space
    $header =~ s/\|\|/\|/g; # Replace || with |
    $header =~ s/^pdb\|//;  # Remove pdb|
    $header =~ s/^sp\|//;   # Remove sp|
    $header =~ s/^pir\|//;  # Remove pir|
    $header =~ s/^prf\|//;  # Remove prf|
    $header =~ s/\|/_/g;    # Replace | with _
    
    $header = ">tcr${header}"; # Put >tcr on the start
    return($header);
}

sub PrintFaa
{
    my($out, $header, $sequence, $minlen) = @_;

    if((length($sequence) >= $minlen) &&
       (($header =~ /tcr/i) ||
        ($header =~ /t-cell\s+receptor/i) ||
        ($header =~ /t-cell-receptor/i) ||
        ($header =~ /t\s+cell\s+receptor/i)) &&
       !($header =~ /antibod/i) &&
       !($header =~ /hybrid/i) &&
       !($header =~ /hypothetical/i))
    {
        $header = FixHeader($header);
        
        print $out "$header\n";
        while(length($sequence))
        {
            my $this = substr($sequence, 0, 80);
            print $out "$this\n";
            $sequence = substr($sequence, 80);
        }
    }
    else
    {
        printf STDERR "Info: Rejected $header%s\n",
            (length($sequence < $minlen)?" (length)":"");
    }
}
