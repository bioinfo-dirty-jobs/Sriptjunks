#!/bin/bash
input="accepted_hits.bam"
lista_nomi=(`ls -d */`)
tm=$(date);
echo "# $tm";
for i in "${lista_nomi[@]}";do Nam=$(echo $i|sed "s/_ens\///g");
echo samtools sort  $i$input   $i$Nam"_sort";done


