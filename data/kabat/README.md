# Obtain Kabat sequence files containing TCRs
kabatDir="/data/kabat/fixlen/unpacked/ALL"
(cd $kabatDir; grep 'TCR BETA CHAIN VARIABLE REGION'  * | grep -v PSEUDO) > beta.txt
(cd $kabatDir; grep 'TCR ALPHA CHAIN VARIABLE REGION' * | grep -v PSEUDO) > alpha.txt 
(cd $kabatDir; grep 'TCR GAMMA CHAIN VARIABLE REGION' * | grep -v PSEUDO) > gamma.txt 
(cd $kabatDir; grep 'TCR DELTA CHAIN VARIABLE REGION' * | grep -v PSEUDO) > delta.txt 

# Now obtain the sequences in FASTA format:
for type in alpha beta gamma delta
do
   ./getkabatseqs.pl $type.txt > $type.faa
done

# ***********************************************************************
# Now align using clustal-omega and obtain the alignments in FASTA format
# with the extension .aln
# ***********************************************************************

# Now remove the breaks to look at the output
for type in alpha beta gamma delta
do
   ./unbreak.pl $type.aln > $type.alignment
done

# ***********************************************************************
# Based on the .alignment files, select a threshold for the number of
# missing AAs that we will accept at each end
# ***********************************************************************
alphaN=1
alphaC=3
betaN=3
betaC=4
gammaN=8
gammaC=4
deltaN=3
deltaC=7

# Keep the sequences that are OK with those boundaries
./keepseqs.pl alpha $alphaN $alphaC alpha.alignment > alpha_final.faa
./keepseqs.pl beta  $betaN  $betaC  beta.alignment  > beta_final.faa
./keepseqs.pl gamma $gammaN $gammaC gamma.alignment > gamma_final.faa
./keepseqs.pl delta $deltaN $deltaC delta.alignment > delta_final.faa

# combine into a single FAA file for use with Fasta
cat *_final.faa >all.faa

#alias fasta=/usr/local/apps/fasta33/fasta33_t
#fasta -q -p -b 10 -d 0 ../tcrseqs.faa all.faa 2 

# Run each of our old TCR sequences against all.faa to see which type categories they go into
./RunFasta.pl ../tcrseqs.faa all.faa alpha_old.faa beta_old.faa gamma_old.faa delta_old.faa

# Now add those to the Kabat sequences
for type in alpha beta gamma delta
do
   cat ${type}_final.faa ${type}_old.faa > ${type}_newold.faa
done

# ***********************************************************************
# Now align the *_newold.faa files using clustal-omega and obtain the
# alignments in FASTA format with the extension .aln
# ***********************************************************************


# Now remove the breaks to look at the output
for type in alpha beta gamma delta
do
   ./unbreak.pl ${type}_newold.aln > ${type}_newold.alignment
done

# ***********************************************************************
# Based on the .alignment files, find the points at which the old (Kabat)
# sequences start and finish
# emacs alpha_newold.alignment alpha.alignment 
# ***********************************************************************
#023014
alphaN=508
alphaC=866
# 028826
betaN=361
betaC=791
# 011027
gammaN=53
gammaC=241
# 035188
deltaN=168
deltaC=390


# Trim the sequences to the boundaries
./trimseqs.pl $alphaN $alphaC alpha_newold.alignment > alpha_newold_trimmed.alignment
./trimseqs.pl $betaN  $betaC  beta_newold.alignment  > beta_newold_trimmed.alignment
./trimseqs.pl $gammaN $gammaC gamma_newold.alignment > gamma_newold_trimmed.alignment
./trimseqs.pl $deltaN $deltaC delta_newold.alignment > delta_newold_trimmed.alignment

# Remove any sequences that don't give the coverage we used earlier
alphaN=1
alphaC=3
betaN=3
betaC=4
gammaN=9
gammaC=4
deltaN=3
deltaC=7

# Keep the sequences that are OK with those boundaries
./keepseqs.pl alpha $alphaN $alphaC alpha_newold_trimmed.alignment > alpha_newold_final.faa
./keepseqs.pl beta  $betaN  $betaC  beta_newold_trimmed.alignment  > beta_newold_final.faa
./keepseqs.pl gamma $gammaN $gammaC gamma_newold_trimmed.alignment > gamma_newold_final.faa
./keepseqs.pl delta $deltaN $deltaC delta_newold_trimmed.alignment > delta_newold_final.faa

# And finally concatenate these files
cat *_newold_final.faa > final.faa



