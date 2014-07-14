#!/bin/bash
tm=$(date)
echo "#$tm"

sito=$1;

mitico=(`find  $sito -name "*1.fastq"`)

for i in "${mitico[@]}";do pif=(`echo $i|sed 's/1.fastq/2.fastq/g'`);
nome=(`echo $pif|awk -F "/" '{print $2}'|sed 's/_2.fastq//g'`);
echo "~/software/defuse-0.6.1/scripts/defuse.pl  -c ~/software/defuse-0.6.1/scripts/config.txt -1" $i "-2" $pif " -o Fusion_defuse"$nome >> defuse_tmp.sh;
echo "~/software/tophat-2.0.9.Linux_x86_64/tophat2 -p 14 --segment-length 25 -o Fusion_tophat"$nome   " --fusion-search   -p 14 --keep-fasta-order --bowtie1  --no-coverage-search --mate-std-dev 80 --max-intron-length 100000 --fusion-min-dist 100000 --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM ~/databases/bowtieDatabase/hg19 $i $pif;" >> Tophat_run_fusion.sh;


echo "~/software/ericscript-0.4.3/ericscript.pl  -db ~/databases/ericscriptdb/ --bwa_aln --nthreads 10 -r /home/sbsuser/databases/bwa_index/hg19.fa -o"   "Eric_"$nome  $i $pif >> Fusion_eric.sh;

echo "python2.7  ~/.local/bin/chimerascan_run.py -p 12 --quals solexa ~/database/chimerascan/" $i   $pif  "chimearascan_"$nome >> Chimerascan_run.sh;
echo "fusioncatcher -d /home/master/database/human_fusioncacther2/ -i"  $i","$pif  "-o  Fusioncatcher_"$nome >> Fusioncatcher_run.sh;


done
