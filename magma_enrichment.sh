#!/bin/bash

myDir=$(pwd)                              #Locate current working directory
magma=($(${myDir}/read_param.sh magma))   #read the location of the magma exe
clozuk=($(${myDir}/read_param.sh clozuk)) #read the location for the SCZ2+CLOZUK GWAS results
asd=($(${myDir}/read_param.sh asd))
depressionCTG=($(${myDir}/read_param.sh depressionCTG))
intelligenceCTG=($(${myDir}/read_param.sh intelligenceCTG))
neuroticismCTG=($(${myDir}/read_param.sh neuroticismCTG))
igap=($(${myDir}/read_param.sh igap))
BP=($(${myDir}/read_param.sh BP))
HD=($(${myDir}/read_param.sh HD))

#for dir in gene_sets; do
#mkdir -p ${myDir}/../output/${dir}
    #list files in directory
# Corrected line to list files in the directory (removed double quotes)

#change the directory for the files 1) baseline 2) 20.1 3)103.2
files=(${myDir}/../gene_sets/103.2/*)

# Loop through files to do the analysis
for file in "${files[@]}"; do
    for gwas in $clozuk $asd $depressionCTG $intelligenceCTG $neuroticismCTG  $igap $BP $HD; do
        if [ -f "$file" ]; then
        #change the output directory for 1)baseline 2)20.1 3)103.2
            output=${myDir}/../output/103.2/
            g="$(basename -- $gwas)"
            f="$(basename -- $file)"
            out=$output$g$f
            ${magma} --gene-results ${gwas} --set-annot "$file" col=2,1 --out ${out}
            fi
    done
done

