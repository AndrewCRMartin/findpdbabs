findpdbabs
==========

A new program to find PDB files containing antibodies

Prerequisites
-------------

- A local mirror of the PDB (you can use our `ftpmirror` script for this)
- BiopTools installed and available in the path
- Legacy blast installed and in the path

Configuration
-------------

Edit the file `findpdbabs.conf` to specify:

- `pdbdir` - the location of your PDB mirror
- `faadir` - the location of the Fasta equivalent of the PDB (the data
  will be created during the processing - you do not need to have
  pre-created it)
- `abpdbdir` - the location of a directory containing confirmed known
  PDB files of antibodies
- `dataroot` - root directory of where you wish to store Blast database
  file, database files etc.

Running the software
--------------------

### First run only

**If you are running the software for the first time**, you need to
prepare the data files for `findpdbabs` to use.

You will need

- A directory containing known antibody structure Fvs (e.g. from AbDb)


Simply ensure you are
ine the `dataprep/findpdbabs` directory and type:

```
   ./builddata.sh
```

This obtains non-redundant FASTA files of the known antibodies,
non-redundant TCR sequences and a file containing them both.

### Subsequent runs

**Once you have done that (and for subsequent runs**, simply type:

```
   ./update.sh
```

This first creates a sequence database based on PDB files. It uses
`getpdbseqs.pl` to creates a sequence file in `faadir` for each PDB
file not already processed.

It then runs runs the program (`./findpdbabs.pl`) to identify new PDB
that contain an antibody.

Algorithm
---------

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

