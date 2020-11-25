pdbdir=[%pdbdir%]
faadir=[%faadir%]

# Get rid of all built-in SUFFIX rules
.SUFFIXES : 

# Delete targets which can't be properly built
# Hmmmm this doesn't actually work cos the target is not placed in the
# current directory
.DELETE_ON_ERROR :

# Continue after all errors
.IGNORE :

# Source and destination directories
vpath %.ent $(pdbdir)
vpath %.faa $(faadir)

# Get a directory listing of all PDB files
sources := $(wildcard $(pdbdir)/*.ent)
# Strip the path
stems1  := $(notdir $(sources))
# Strip the extension
stems   := $(basename $(stems1))
# Add the .faa extension to create list of faa targets
faatargets := $(addsuffix .faa, $(stems))

# Phony rule to build all the targets
all : faa

faa : $(faatargets)

%.faa : %.ent
	(cd $(faadir); pdb2pir -s -c -f $< > $@)

