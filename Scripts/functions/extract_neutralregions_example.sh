#!/bin/bash
# Extract sequences of short introns at position 8-30

while read line ; do bedtools getfasta -fi ${line}_Chr2L_wh.fa -bed /Volumes/Burcin/Reference/Introns/final_SI/bed_files/uniq/neutral_regions/SI_2L_complement_neutral.bed > last_SI/${line}_Chr2L_SI_complement.fa ; done < /Volumes/Burcin/Drosophila_data/Dmel/individuals

while read line ; do bedtools getfasta -fi ${line}_Chr2L_wh.fa -bed /Volumes/Burcin/Reference/Introns/final_SI/bed_files/uniq/neutral_regions/SI_2L_noncomplement_neutral.bed > last_SI/${line}_Chr2L_SI_noncomplement.fa ; done < /Volumes/Burcin/Drosophila_data/Dmel/individuals
