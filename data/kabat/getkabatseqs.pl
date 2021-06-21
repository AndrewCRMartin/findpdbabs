#!/usr/bin/perl
use strict;

my $kabatDir = "/data/kabat/fixlen/unpacked/ALL";

while(<>)
{
    my($id) = split(/:/);
    if($id =~ /^\d/)
    {
        my($name, $species, $definition, $sequence) = ReadKabat($id, $kabatDir);
        if(length($sequence) > 70)
        {
            my $label = ">$id|$name|$species|$definition";
            $label =~ s/\s/_/g;
            print "$label\n";
            print "$sequence\n";
        }
    }
}

sub ReadKabat
{
    my($code, $kabatDir) = @_;

    my $name       = '';
    my $species    = '';
    my $definition = '';
    my $sequence   = '';
    
    my $file = "$kabatDir/$code";
    if(open(my $fp, '<', $file))
    {
        while(<$fp>)
        {
            chomp;
            if(/^DEFINI\s+(.*)$/)
            {
                $definition = $1;
            }
            elsif(/^AANAME\s+(.*)$/)
            {
                $name = $1;
            }
            elsif(/^SPECIE\s+(.*)$/)
            {
                $species = $1;
            }
            elsif(/^SEQTPA/)
            {
                my(@fields) = split();
                if(defined($fields[5]))
                {
                    if($fields[5] ne '-')
                    {
                        $sequence .= $fields[5];
                    }
                }
            }
        }
        close($fp);
    }
    else
    {
        printf STDERR "Can't read $file\n";
    }

    return($name, $species, $definition, $sequence);
}
