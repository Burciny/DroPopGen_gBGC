#!/bin/bash

awk -F "A" '{print NF-1}' $1 > A
awk -F "T" '{print NF-1}' $1 > T
awk -F "G" '{print NF-1}' $1 > G
awk -F "C" '{print NF-1}' $1 > C
awk -F "N" '{print NF-1}' $1 > N
paste A T G C N > $2
rm A T G C N


