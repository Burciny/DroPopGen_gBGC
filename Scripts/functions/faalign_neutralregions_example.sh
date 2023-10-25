#!/bin/bash

# Change header with individual names

while read line ; do sed "s/^\(>\).*$/\1${line}/1" ${line}_Chr2L_SI_complement.fa > ${line}_Chr2L_SI_complement_hc.fa ; done < /Volumes/Burcin/Drosophila_data/Dmel/individuals
while read line ; do sed "s/^\(>\).*$/\1${line}/1" ${line}_Chr2L_SI_noncomplement.fa > ${line}_Chr2L_SI_noncomplement_hc.fa ; done < /Volumes/Burcin/Drosophila_data/Dmel/individuals

# Combine same introns from different individuals
a=1
b=1
while [ $a -le 4699 ] ; do 
	while read line ; do
		sed -n ${a},$((a+1))p Chr2L/fasta/last_SI/${line}_Chr2L_SI_complement_hc.fa >> Alignments_SI_final/neutral/2L/SI_2L_complement_NR_${b}.fa
	done < /Volumes/Burcin/Drosophila_data/Dmel/individuals
	a=$((a+2))
	b=$((b+1))
done


a=1
b=1
while [ $a -le 4583 ] ; do 
	while read line ; do
		sed -n ${a},$((a+1))p Chr2L/fasta/last_SI/${line}_Chr2L_SI_noncomplement_hc.fa >> Alignments_SI_final/neutral/2L/SI_2L_noncomplement_NR_${b}.fa
	done < /Volumes/Burcin/Drosophila_data/Dmel/individuals
	a=$((a+2))
	b=$((b+1))
done
