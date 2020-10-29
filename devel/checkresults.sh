#!/bin/bash

pdbs=`awk -F_ '{print $1}' abs.out | sort -u`

for pdb in $pdbs
do
    file="/serv/data/pdb/pdb${pdb}.ent"
    match=`egrep 'COMPND.*MOLECULE' $file`

    if [[ $match =~ "ANTIBODY" ]] ||
           [[ $match =~ " FV " ]] ||
           [[ $match =~ " FV;" ]] ||
           [[ $match =~ " FV," ]] ||
           [[ $match =~ " VH;" ]] ||
           [[ $match =~ " VH-" ]] ||
           [[ $match =~ " VL;" ]] ||
           [[ $match =~ " VL-" ]] ||
           [[ $match =~ " MAB " ]] ||
           [[ $match =~ " MAB;" ]] ||
           [[ $match =~ " IGM " ]] ||
           [[ $match =~ " IGG " ]] ||
           [[ $match =~ ",IGG " ]] ||
           [[ $match =~ " IGG;" ]] ||
           [[ $match =~ "IGG1" ]] ||
           [[ $match =~ "IGG2A" ]] ||
           [[ $match =~ "IGA1" ]] ||
           [[ $match =~ " SCFV" ]] ||
           [[ $match =~ " VHH" ]] ||
           [[ $match =~ "ANTI-" ]] ||
           [[ $match =~ "DIABODY" ]] ||
           [[ $match =~ "NANOBODY" ]] ||
           [[ $match =~ "DSFV" ]] ||
           [[ $match =~ "SINGLE-CHAIN VARIABLE" ]] ||
           [[ $match =~ "SINGLE CHAIN VARIABLE" ]] ||
           [[ $match =~ "IMMUNOGLOBULIN E " ]] ||
           [[ $match =~ "IMMUNOGLOBULIN E;" ]] ||
           [[ $match =~ "IMMUNOGLOBULIN A " ]] ||
           [[ $match =~ "IMMUNOGLOBULIN A;" ]] ||
           [[ $match =~ "IMMUNOGLOBULIN G " ]] ||
           [[ $match =~ "IMMUNOGLOBULIN G;" ]] ||
           [[ $match =~ "IMMUNOGLOBULIN D " ]] ||
           [[ $match =~ "IMMUNOGLOBULIN D;" ]] ||
           [[ $match =~ "IMMUNOGLOBULIN E " ]] ||
           [[ $match =~ "IMMUNOGLOBULIN E;" ]] ||
           [[ $match =~ "IMMUNOGLOBULIN M " ]] ||
           [[ $match =~ "IMMUNOGLOBULIN M;" ]] ||
           [[ $match =~ "IMMUNOGLOBULIN H " ]] ||
           [[ $match =~ "IMMUNOGLOBULIN H;" ]] ||
           [[ $match =~ "IMMUNOGLOBULIN L " ]] ||
           [[ $match =~ "IMMUNOGLOBULIN L;" ]] ||
           [[ $match =~ "IMMUNOGLOBULIN LAMBDA" ]] ||
           [[ $match =~ "IMMUNOGLOBULIN HEAVY" ]] ||
           [[ $match =~ "IMMUNOGLOBULIN KAPPA" ]] ||
           [[ $match =~ "VH DOMAIN" ]] ||
           [[ $match =~ "BENCE JONES" ]] ||
           [[ $match =~ "BENCE-JONES" ]] ||
           [[ $match =~ "HEAVY CHAIN" ]] ||
           [[ $match =~ "HEAVY-CHAIN" ]] ||
           [[ $match =~ "LIGHT CHAIN" ]] ||
           [[ $match =~ "(L CHAIN" ]] ||
           [[ $match =~ " L CHAIN" ]] ||
           [[ $match =~ "(H CHAIN" ]] ||
           [[ $match =~ " H CHAIN" ]] ||
           [[ $match =~ "LAMBDA DIMER" ]] ||
           [[ $match =~ "LAMBDA DOMAIN" ]] ||
           [[ $match =~ "LAMBDA CHAIN" ]] ||
           [[ $match =~ "KAPPA CHAIN" ]] ||
           [[ $match =~ "KAPPA DOMAIN" ]] ||
           [[ $match =~ "LIGHT CHAIN" ]] ||
           [[ $match =~ "HEAVY CHAIN" ]] ||
           [[ $match =~ "HEAVY DOMAIN" ]] ||
           [[ $match =~ "KAPPA-4" ]] ||
           [[ $match =~ "IG KAPPA" ]] ||
           [[ $match =~ "IG LAMBDA" ]] ||
           [[ $match =~ "IG GAMMA" ]] ||
           [[ $match =~ " FAB" ]]; then
#        echo "$pdb is an ANTIBODY"
        echo -n "."
    else
        echo ""
        echo $pdb
        echo $match
        echo ""
        echo "--------------------------------------"
        echo ""
    fi
done
echo ""

echo -n "Total number of potential antibodies: "
echo $pdbs | wc -w
