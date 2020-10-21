#!/usr/bin/perl
use strict;

my $pdbdir="/serv/data/pdb";
my $allseq="alltcr.faa";
my $repseq="tcrseqs.faa";
my $tcrlist="tcrpdbs.dat";
$::pdbcount="pdbcount";

unlink $allseq;

my @pdbfiles = ReadFile($tcrlist);
@pdbfiles = CheckFiles($pdbdir, @pdbfiles);

print STDERR scalar(@pdbfiles) . "\n";

print STDERR "Concatenating files";
foreach my $file (@pdbfiles)
{
    my $pdbcode = $file;
    $pdbcode =~ s/.*\///;
    $pdbcode =~ s/\..*?$//;
    $pdbcode =~ s/^pdb/tcr/;
    `pdb2pir -c -f -l ${pdbcode}_ ${pdbdir}/$file >>$allseq`;
     print STDERR '.';
}
print STDERR "done\n";
# Grab and compile CD-HIT if not there
if( ! -d "cdhit" )
{
    print STDERR "Building cd-hit\n";
    `git clone git\@github.com:weizhongli/cdhit.git`;
    `(cd cdhit; make)`;
    print STDERR "done\n";
}

print STDERR "Running CD-HIT...";
`cdhit/cd-hit -c 0.6 -n 3 -i $allseq -o $repseq`;
print STDERR "done\n";

#SplitFastaFiles($repseq, $tplDir);


#*************************************************************************
sub CheckFiles
{
    my($pdbdir, @pdbfiles) = @_;
    my @newfiles = ();
    print STDERR "Checking for files with exactly 2 full-length chains";
    foreach my $pdbfile (@pdbfiles)
    {
        my $filename = "$pdbdir/$pdbfile";
        if( -e $filename)
        {
            print STDERR '.';
            my $exe = "$::pdbcount -c $filename";
            chomp $exe;
            my $result = `$exe`;
            my @lines = split(/\n/, $result);
            my $ok = 1;
            if(scalar(@lines) == 3)  # we have 2 chains plus totals line
            {
                for(my $lineCount=0; $lineCount<2; $lineCount++)
                {
                    my @fields = split(/\s+/, $lines[$lineCount]);
                    if($fields[3] < 100)
                    {
                        $ok = 0;
                    }
                }
                if($ok)
                {
                    push @newfiles, $pdbfile;
                    print "$pdbfile\n" if(defined($::v));
                }
            }
        }
    }
    print STDERR "done\n";

    return(@newfiles);
}

#*************************************************************************
sub ReadFile
{
    my($file) = @_;
    my @files = ();
    if(open(my $fp, '<', $file))
    {
        while(<$fp>)
        {
            chomp;
            $_ = "\L$_"; # downcase
            $_ = "pdb$_" . ".ent";
            push @files, $_;
        }
        close $fp;
    }
    return(@files);
}

