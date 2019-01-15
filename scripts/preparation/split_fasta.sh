#!/bin/bash 
csplit -s -z $1 '/>/' '{*}'
for i in xx* ; do \
    n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
    mv "$i" "$n.fa" ; \
done
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
    if [ -f $chrom.fa ]
    then
        mv $chrom.fa chr$chrom.fa
    fi
done
