#!/bin/sh

echo "Indel coding log" > multi-indel-code.log.txt
echo >> multi-indel-code.log.txt

while [ $# -gt 0 ]
do

	echo "locus -- " $1
	./code-individual-introns.pl 2500orthologs-introns-filtered $1 jarvis-taxa-authority.txt indel-data-2516orthologs-introns-filtered
	cat indel-data-2516orthologs-introns-filtered/intronic_locus_$1_indels.log.txt >> multi-indel-code.log.txt
	
	shift
	
done
