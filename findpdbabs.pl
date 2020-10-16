#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    
#   File:       
#   
#   Version:    
#   Date:       
#   Function:   
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin, UCL, 2011
#   Author:     Dr. Andrew C. R. Martin
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
# Or if we have a bin directory and a lib directory
#use Cwd qw(abs_path);
#use FindBin;
#use lib abs_path("$FindBin::Bin/../lib");

use strict;
use config;

my $configFile = 'findpdbabs.conf';

#*************************************************************************
my %config = config::ReadConfig($configFile);

$::np        = 2 if(!defined($::np));
$::minLen    = 80;
$::minId     = 0.3;
$::maxEValue = 0.01;
$|=1;

my $faaDir    = $config{'faadir'};
my $dbmFile   = $config{'dbmfile'};
my $seqFile   = $config{'seqfile'};
my $tplDir    = $config{'tpldir'};
my %processed = ();

my $nSeqs = BuildSequenceFile($faaDir, $dbmFile, $seqFile);
if($nSeqs > 0)
{
    my $tmpDir = RunBlast($tplDir, $seqFile, $nSeqs);
    if($tmpDir eq '')
    {
        print STDERR "Error: BLAST failed\n";
        exit 1;
    }
    my @abs = FindAbs($tmpDir);
    unlink $tmpDir;
}
else
{
    exit 1;
}

#*************************************************************************
sub FindAbs
{
    my($blastDir)  = @_;

    my %lengths    = ();
    my %evalues    = ();
    my %ids        = ();
    my %chainTypes = ();
    my %positives  = ();
    my @blastFiles = GetFileList($blastDir, '.out');
    
    foreach my $blastFile (@blastFiles)
    {
        print STDERR "Parsing $blastFile...";
        ParseBlast($blastFile, \%lengths, \%evalues, \%ids,
                   \%positives, \%chainTypes);
        print STDERR "done\n";
    }

    my @labels = keys %lengths;
#    ReverseSearch
    OutputAbs(\@labels, \%lengths, \%evalues, \%ids, \%positives, \%chainTypes);
}

sub OutputAbs
{
    my($aLabels, $hLengths, $hEvalues, $hIds, $hPositives, $hChainTypes) = @_;

#    foreach my $label (reverse sort {$$hIds{$a} <=> $$hIds{$b}} @$aLabels)
    foreach my $label (reverse sort {$$hPositives{$a} <=> $$hPositives{$b}} @$aLabels)
    {
        printf("$label: ChainType: $$hChainTypes{$label} E: %6.2g ID: %.2f Pos: %.2f Len: %3d\n",
            $$hEvalues{$label}, $$hIds{$label}, $$hPositives{$label}, $$hLengths{$label});
    }
}

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
    
