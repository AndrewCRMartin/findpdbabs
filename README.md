findpdbabs
==========

A new program to find PDB files containing antibodies

Prerequisites
-------------

- A local mirror of the PDB
- BiopTools installed and available in the path
- Legacy blast installed and in the path
- A directory containing known antibody structure Fvs (e.g. from AbDb)

Configuration
-------------

Edit the file `findpdbabs.conf` to specify:

- `faadir` - the location of the Fasta equivalent of the PDB
- `abpdbdir` - the location of a directory containing confirmed known PDB files of antibodies

Edit the file `Makefile.pdbseq` to specify:

- `faadir` - the location of the Fasta equivalent of the PDB (the same as above)
- `pdbdir` - the location of the local PDB mirror

Running the software
--------------------

This requires various steps:

1. Prepare the data files for findpdbabs to use. Simply type: `builddata.sh`. This obtains non-redundant FASTA files of the known antibodies, non-redundant TCR sequences and a file containing them both.

2. First we need a sequence database based on PDB files. This is done
by the command: `make -f Makefile.pdbseq` which will simply create
the sequence files for any PDB files that don't have them

3. Now run the program to identify new PDB sequence files and identify
them as containing an antibody: `./findpdbabs.pl`

** For future runs, you just repeat steps 2 and 3. This will simply
   identify new antibody sequences **

The algorithm is to identify PDB sequences that haven't previously
been examined (a DBM file is used to keep track of this) and scan each
of the template sequences against these. If the match is over >= 80
residues and had a sequence ID of >= 0.3 with an e-value <= 0.01, then
the hit is retained. When each template is run with BLAST, the results
will replace older results for a given new sequence hit if the new
match is better.

The next stage is a BLAST search of each hit against a database of
antibody and TCR sequences. If the best hit is an antibody, the
sequence is kept, if it is a TCR (or doesn't match anything) it is
rejected.

