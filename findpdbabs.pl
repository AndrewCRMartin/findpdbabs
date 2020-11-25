#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    findpdbabs
#   File:       findpdbabs.pl
#   
#   Version:    V1.0
#   Date:       27.10.20
#   Function:   Find PDB files containing antibodies
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

$::np        = 2 if(!defined($::np));
$::minLen    = 100;
$::minId     = 0.3;
$::maxEValue = 1e-20;
$|=1;

UsageDie() if(defined($::h) ||
              (!defined($::check) && (scalar(@ARGV) != 1)));
$|=1;

my $outFile    = shift(@ARGV);

my $faaDir     = $config{'faadir'};
my $dbmFile    = $config{'dbmfile'};
my $seqFile    = $config{'seqfile'};
#my $tplDir     = "${exeDir}/" . $config{'tplsubdir'};
my $tplDir     = $config{'abtpldir'};
my $tcrAbsFile = $config{'tcrabsfile'};
my %processed  = ();

if(defined($::check))
{
    CheckOneSequence($::check, $seqFile, $tcrAbsFile);
    exit 0;
}

if(open(my $outFp, '>', $outFile))
{
    my $nSeqs = BuildSequenceFile($faaDir, $dbmFile, $seqFile);
    if($nSeqs > 0)
    {
        my $tmpDir = RunBlast($tplDir, $seqFile, $nSeqs);
        if($tmpDir eq '')
        {
            print STDERR "Error: BLAST failed\n";
            exit 1;
        }
        my @abs = FindAndPrintAbs($outFp, $tmpDir, $seqFile, $tcrAbsFile);
        if(defined($::d))
        {
            print "BLAST output is in $tmpDir\n";
        }
        else
        {
            unlink $tmpDir;
        }
    }
}
else
{
    print STDERR "Error: Unable to write to output file ($outFile)\n";
    exit 1;
}


#*************************************************************************
sub FindAndPrintAbs
{
    my($outFp, $blastDir, $seqFile, $tcrAbsFile)  = @_;

    my %lengths    = ();
    my %evalues    = ();
    my %ids        = ();
    my %chainTypes = ();
    my %positives  = ();
    my @blastFiles = GetFileList($blastDir, '.out');
    
    foreach my $blastFile (@blastFiles)
    {
        print STDERR "Parsing $blastFile...";
        ParseBlast("$blastDir/$blastFile", \%lengths, \%evalues, \%ids,
                   \%positives, \%chainTypes);
        print STDERR "done\n";
    }

        
    
    my @labels = keys %lengths;

    printf STDERR "\n\nBLAST searches found %d possible antibodies\n\n",
        scalar(@labels);


    @labels = ReverseSearch($tcrAbsFile, $seqFile, @labels);
    OutputAbs($outFp, \@labels, \%lengths, \%evalues, \%ids,
              \%positives, \%chainTypes);
}


#*************************************************************************
sub ReverseSearch
{
    my($tcrAbsFile, $seqFile, @labels) = @_;
    my $logFp;
    my $nHits  = scalar(@labels);
    my $hitNum = 1;
    
    if(defined($::l) && ($::l ne ''))
    {
        if(!open($logFp, '>', $::l))
        {
            printf STDERR "Warning: Unable to open log file ($::l) - log going to standard error\n";
            $logFp = 0;
        }
    }
    
    my @newlabels = ();
    if(BuildBlastDB($tcrAbsFile))
    {
        foreach my $label (@labels)
        {
            print STDERR "($hitNum/$nHits): ";
            if(BlastCheck($label, $seqFile, $tcrAbsFile))
            {
                push @newlabels, $label;
                PrintLog($logFp, "Retained Ab: $label");
            }
            else
            {
                PrintLog($logFp, "Rejected TCR: $label");
            }
            $hitNum++;
        }
    }

    close($logFp) if($logFp);
    
    return(@newlabels);
}


#*************************************************************************
sub PrintLog
{
    my($logFp, $msg) = @_;
    if(defined($::v) || defined($::l))
    {
        if($logFp)
        {
            print $logFp "$msg\n";
        }
        else
        {
            print STDERR "$msg\n";
        }
    }
}


#*************************************************************************
sub BlastCheck
{
    my($label, $seqFile, $tcrAbsFile) = @_;
    my $isAntibody = 0;

    print STDERR "Reverse BLAST for label ${label}...";
    my $testFile = GetSequence($seqFile, $label);
    if($testFile ne '')
    {
        my $outFile = $testFile . $label . ".out";
        my $exe = "blastall -F F -e $::maxEValue -z 100000 -P 1 -p blastp -a $::np -d $tcrAbsFile -i $testFile -o $outFile";
        `$exe`;
        unlink($testFile);
        $isAntibody = CheckAntibody($outFile);
        if(defined($::d))
        {
            if($::d > 1)
            {
                print STDERR "\n   Result: $outFile\n";
            }
        }
        else
        {
            unlink($outFile);
        }
    }
    print STDERR "done\n";
    
    return($isAntibody);
}


#*************************************************************************
sub CheckAntibody
{
    my($blastFile) = @_;
    my $retval = 0;
    if(open(my $fp, '<', $blastFile))
    {
        my $inData = 0;
        while(<$fp>)
        {
            if(/^Sequences producing/)
            {
                $_ = <$fp>;  # Get blank line
                $_ = <$fp>;  # Get first hit
                if(/^tcr/)
                {
                    $retval = 0; # It's a TCR
                }
                elsif(/^[0-9]/)
                {
                    $retval = 1; # It's a PDB identifier (i.e. and antibody)
                }
                elsif(/No hits/)
                {
                    $retval = 0; # It's not found so not an antibody
                }
                else
                {
                    $retval = 0; # We don't know what it is!
                }
                last;
            }
        }
        close($fp);
    }
    return($retval);
}


#*************************************************************************
sub GetSequence
{
    my($seqFile, $label) = @_;
    my $tmpFile = '';
    
    if(open(my $fp, '<', $seqFile))
    {
        my $inData = 0;
        my $result = '';
        while(<$fp>)
        {
            if(/^\>$label/)
            {
                $result .= $_;
                $inData = 1;
            }
            elsif(/^\>/)
            {
                last if($inData);
                $inData = 0;
            }
            elsif($inData)
            {
                $result .= $_;
            }
        }
        close($fp);

        if(length($result))
        {
            $tmpFile = "/var/tmp/fpatest_$$" . "__" . time();
            if(open(my $out, '>', $tmpFile))
            {
                print $out $result;
                close $out;
            }
            else
            {
                $tmpFile = '';
            }
        }
    }
    return($tmpFile);
}


#*************************************************************************
sub OutputAbs
{
    my($outFp, $aLabels, $hLengths, $hEvalues, $hIds, $hPositives, $hChainTypes) = @_;

#    foreach my $label (reverse sort {$$hIds{$a}*1.0 <=> $$hIds{$b}*1.0} @$aLabels)
    foreach my $label (reverse sort {$$hPositives{$a}*1.0 <=> $$hPositives{$b}*1.0} @$aLabels)
    {
        printf($outFp "$label: ChainType: $$hChainTypes{$label} E: %6.2g ID: %.2f Pos: %.2f Len: %3d\n",
            $$hEvalues{$label}, $$hIds{$label}, $$hPositives{$label}, $$hLengths{$label});
    }
}


#*************************************************************************
sub ParseBlast
{
    my($blastFile, $hLengths, $hEvalues, $hIds, $hPositives, $hChainTypes) = @_;
    if(open(my $fp, '<', $blastFile))
    {
        my $id  = 0;
        my $pos = 0;
        my $len = 0;
        my $e   = 0;
        my $label = '';
        my $chainType = $blastFile;
        $chainType =~ s/.*_//;
        $chainType =~ s/\..*//;
        
        while(<$fp>)
        {
            chomp;
            if(/\>(.*)/)
            {
                $label = $1;
            }
            elsif(/Expect\s+=\s+(.*?)\,/)
            {
                $e = $1;
            }
            elsif(/Identities\s+=\s+(\d+)\/(\d+).*Positives\s+=\s+(\d+)\/.*/)
            {   #                    id     len                    pos
                $len = $2;
                $id  = $1/$len;
                $pos = $3/$len;
            }
            elsif(/Query/)
            {
                if($label ne '')
                {
                    # print "$label: ChainType: $chainType E: $e ID: $id Pos: $pos Len: $len\n";
                    UpdateStats($label, $hLengths, $hEvalues, $hIds, $hPositives, $hChainTypes,
                                        $len,      $e,        $id,   $pos,        $chainType);
                    $id  = 0;
                    $pos = 0;
                    $len = 0;
                    $e   = 0;
                    $label = '';
                }
            }
        }
        close $fp;
    }
    else
    {
        print STDERR "Error: Unable to read $blastFile\n";
    }
}


#*************************************************************************
sub UpdateStats
{
    my($label, $hLengths, $hEvalues, $hIds, $hPositives, $hChainTypes,
               $len,      $e,        $id,   $pos,        $chainType) = @_;

    $len = int($len);
    
    if(($len > $::minLen) &&
       ($id  > $::minId))
    {
        {
            if(!defined($$hLengths{$label}))
            {
                $$hLengths{$label}    = $len;
                $$hEvalues{$label}    = $e   * 1.0;
                $$hIds{$label}        = $id  * 1.0;
                $$hPositives{$label}  = $pos * 1.0;
                $$hChainTypes{$label} = $chainType;
            }
            else
            {
                if(($e   <  $$hEvalues{$label}) &&
                   ($len >= $$hLengths{$label} - 10))
                {
                    $$hLengths{$label}    = $len;
                    $$hEvalues{$label}    = $e   * 1.0;
                    $$hIds{$label}        = $id  * 1.0;
                    $$hPositives{$label}  = $pos * 1.0;
                    $$hChainTypes{$label} = $chainType;
                }
            }
        }
    }
}


#*************************************************************************
sub RunBlast
{
    my($tplDir, $seqFile, $nSeqs) = @_;

    my $tmpDir = "/var/tmp/fpa_$$" . "_" . time();
    `mkdir $tmpDir` if(! -d $tmpDir);
    
    if(BuildBlastDB($seqFile))
    {
        my @tplFiles = GetFileList($tplDir, '.faa');
        if(defined($::d) && ($::d > 1))
        {
            $#tplFiles = 2;
        }
        
        if(scalar(@tplFiles))
        {
            foreach my $tplFile (@tplFiles)
            {
                print STDERR "Running BLAST with template ${tplFile}...";
                my $exe = "blastall -F F -e $::maxEValue -z 100000 -P 1 -b $nSeqs -p blastp -a $::np -d $seqFile -i $tplDir/$tplFile -o ${tmpDir}/${tplFile}.out";
                `$exe`;
                print STDERR "done\n";
            }
            return($tmpDir);
        }
        print STDERR "No template files found\n";
        return('');
    }
    return(0);
}


#*************************************************************************
sub BuildBlastDB
{
    my ($seqFile) = @_;
    print STDERR "Building BLAST database..." if(defined($::v));

    my $exe = "formatdb -i $seqFile -p T -o F";
    `$exe`;

    if( -f "${seqFile}.pin" )
    {
        print STDERR "done\n" if(defined($::v));
        return(1);
    }
    print STDERR "failed\n" if(defined($::v));
    return(0);
}


#*************************************************************************
sub BuildSequenceFile
{
    my($faaDir, $dbmFile, $seqFile) = @_;
    dbmopen %processed, $dbmFile, 0666 || die "Can't dbopen $dbmFile";

    my $nSeqs = 0;

    print STDERR "Finding new sequence files";
    
    my @fastaFiles = GetFileList($faaDir, '.faa');

    if(defined($::d))
    {
        $#fastaFiles = 1000;
    }
    
    if(scalar(@fastaFiles))
    {
        if(open(my $seqFp, '>', $seqFile))
        {
            foreach my $file (@fastaFiles)
            {
                if(!defined($processed{$file}))
                {
                    my $nseq = AddToSeqFile("$faaDir/$file", $seqFp);
                    if($nseq)
                    {
                        $nSeqs+=$nseq;
                        $processed{$file} = 1;
                        if(defined($::v))
                        {
                            print STDERR "Added $file\n"
                        }
                        else
                        {
                            print STDERR '.';
                        }
                    }
                    else
                    {
                        print STDERR "\nWarning: No sequences added from $file\n";
                    }
                }
            }
            close($seqFp);
            print STDERR "done\n";
            
            if($nSeqs)
            {
                print STDERR "\nAdded $nSeqs new sequences\n";
            }
            else
            {
                print STDERR "No new sequence files\n";
            }
        }
        else
        {
            print STDERR "Error: Unable to open $seqFile for writing\n";
            dbmclose %processed;
            return(0);
        }
    }
    else
    {
        print STDERR "No sequence files found\n";
    }
    
    dbmclose %processed;
    return($nSeqs);
}


#*************************************************************************
sub AddToSeqFile
{
    my($file, $seqFp) = @_;
    my $pdbCode = $file;
    my $nseqs = 0;
    $pdbCode =~ s/.*\///;
    $pdbCode =~ s/^pdb//;
    $pdbCode =~ s/\.faa//;
    my $entry = '';

    if(open(my $in, '<', $file))
    {
        while(<$in>)
        {
            chomp;
            if(/\>Chain(.*)/)
            {
                my $chain = $1;
                my $pdbc = "${pdbCode}_${chain}";
                print $seqFp ">$pdbc\n";
                $nseqs++;
            }
            else
            {
                print $seqFp "$_\n";
            }
        }
        close $in;
        return($nseqs);
    }
    else
    {
        printf STDERR "Warning: Unable to read file $file\n";
    }
    return(0);
}


#*************************************************************************
sub GetFileList
{
    my($dir, $ext) = @_;
    my @files = ();
    if(opendir(my $dirFp, $dir))
    {
        @files = grep(/$ext/, readdir($dirFp));
        closedir($dirFp);
    }
    return(@files);
}


#*************************************************************************
sub CheckOneSequence
{
    my($label, $seqFile, $tcrAbsFile)  = @_;
    if((! -f $seqFile) || (! -f $tcrAbsFile))
    {
        print STDERR "Error: Sequence file ($seqFile) or TCR/Abs file ($tcrAbsFile) doesn't exist\n";
        exit 1;
    }
    $::d = 2; # turn on debug mode
    my $isAb = BlastCheck($label, $seqFile, $tcrAbsFile);
}


#*************************************************************************
sub UsageDie
{
    print <<__EOF;

findpdbabs.pl V1.0  (c) 2020, UCL, Prof. Andrew C.R. Martin

Usage: findpdbabs.pl [-h][-np=n][-d[=n]][-v][-l=logfile] abs.out
-or-   findpdbabs.pl -check=pppp_c
       -h     This help
       -np    Specify number of CPU threads for BLAST [$::np]
       -d     Debug - limits number of sequences used to 1000
                      -d=2 Also limit number of templates for BLAST
                      searches to 2
       -v     Verbose
       -l     Specify log file for listing sequences rejected as TCRs
              rather than sending them to STDERR
       -check Specify a pdbcode_chain (pppp_c) - e.g. 1yqv_L - to
              run just the reverse blast on

findpbdabs.pl uses a pre-calculated set of FASTA sequence files
derived from PDB files to identify antibody chains. It uses a forward
BLAST search using known antibody sequences and a reverse search on
the hits against a set of known antibody and TCR sequuences, rejecting
those that hit a TCR as the top hit.

__EOF

    exit 0;
}
