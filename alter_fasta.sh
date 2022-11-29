#!/bin/bash

organism_name="$1"

if [[ ${organism_name} = "" ]]; then
	echo "no organism name was given"
else
	cat proteoom/*.pep.all.fa | tr -d "\n" | sed 's/>/\n>/g' > ${organism_name}.ensembl.fa 
	cat proteoom/*.fasta | tr -d "\n" | sed 's/>/\n>/g' > ${organism_name}.uniprot.fa 
fi
