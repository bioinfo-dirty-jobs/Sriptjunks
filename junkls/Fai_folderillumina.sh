#!/bin/bash
#script Analisi Dati E allineamento
tm=$(date)
echo "#$tm"
nomi="$1"; #nome del Mese Esperimento e data
sito=$2;


print_usage() {
		echo "This Script genearte all the Scripts we need for Analazed a single run"
		echo "echo ./Execute.all.sh [-name folder]"
		echo "Generate a folder tree for analazyng data"
		echo "Generate the Script_allignment.sh for run Tophat allignment"
		echo "Generate the Fusion_eric.sh for run Ericscript"
		echo "Generate the Fuscatch_run.sh for run FusionCatcher"
		echo "Generate the Tophat_run.sh for run Tophat fusion"
	}

	[ $# != 2 ] && print_usage && exit 1
	


if [ -d /$nomi] ]
then
    echo the directory exists;
    echo $nomi

else
	echo Create the folder
	mkdir "$nomi"
	cd "$nomi"
		mkdir ALIGN
		mkdir Quality
		mkdir POST_ALIGN
		mkdir Fusion_res
		mkdir Fusion_res/Eric
		mkdir Fusion_res/Fucatcher/
		mkdir Fusion_res/Tophatfusion/
		mkdir Fusion_res/chimerascan






fi


#sito="/illumina/runs/130820_SNL129_0048_AD1T4JACXX/Finale/Project_DefaultProject/"

#mitico=(`ls -l  $sito|awk '{print $9}'`)
valr=(`ls -l  $sito|awk '{print $9}'|awk -F "_" '{print $1}'|uniq`);
mitico=(`find  $sito -name "*_R1_*"`)


echo "#!/bin/bash" > Script_allignment.sh;

echo "#PBS -l nodes=0:1:ppn=8,mem=22gb">> Script_allignment.sh;

echo "#$tm">> Script_allignment.sh;
echo "#!/bin/bash" > Tophat_run_fusion.sh;
echo "#PBS -l nodes=0:1:ppn=8,mem=22gb">> Tophat_run_fusion.sh;
echo "#$tm">> Tophat_run_fusion.sh;

echo "#!/bin/bash" > Tophat2_run_hg19.sh;
echo "#PBS -l nodes=0:1:ppn=8,mem=22gb">> Tophat2_run_hg19.sh;
echo "#$tm">> Tophat2_run_hg19.sh;


echo "#!/bin/bash" > Chimerascan_run.sh;

echo "#PBS -l nodes=0:1:ppn=8,mem=22gb">> Chimerascan_run.sh;

echo "#$tm">> Chimerascan_run.sh;









echo "#PBS -l nodes=1:ppn=14,mem=22gb">> Tophat_run_ens.sh;
echo "#$tm">> Tophat_run_ens.sh;
for i in "${mitico[@]}"
do	pet=(`echo $i|sed 's/_R1_/_R2_/g'`);
pref=(`echo $i|awk -F "/" '{print $NF}'|awk -F "_" '{print $1}' `);
echo "~/software/tophat-2.0.7.Linux_x86_64/tophat2 -p 14 --segment-length 25 -o $nomi/ALIGN/$pref -G ~/databases/bowtie2_ens/Homo_sapiens.GRCh37.72.gtf ~/databases/bowtie2_ens/Homo_sapiens.GRCh37.72.dna_sm.toplevel $i  $pet ;"   >> Script_allignment.sh  
echo "~/software/tophat-2.0.9.Linux_x86_64/tophat2  -o $nomi/Fusion_res/Tophatfusion/"$pref"    --fusion-search   -p 14 --keep-fasta-order --bowtie1  --no-coverage-search --mate-std-dev 80 --max-intron-length 100000 --fusion-min-dist 100000 --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM ~/databases/bowtieDatabase/hg19 $i $pet;" >> Tophat_run_fusion.sh 
echo "~/software/tophat-2.0.9.Linux_x86_64/tophat2  -o $nomi/ALIGN/"$pref"   -p15 --segment-length 25 -G ~/databases/tophat/refSeq_hg19.gff   ~/databases/bowtie2Database/hg19  $i  $pet" >> Tophat2_run_hg19.sh 
echo "~/software/tophat-2.0.9.Linux_x86_64/tophat2  -o $nomi/ALIGN/"$pref"   -p15 --segment-length 25 --library-type fr-firststrand -G ~/databases/tophat/refSeq_hg19.gff  ~/databases/bowtie2Database/hg19  $i  $pet" >> Tophat2_run_hg19b.sh 
done



echo "#!/bin/bash" > Fusion_eric.sh;
echo "#!/bin/bash" > Fuscatch_run.sh;
echo "#PBS -l nodes=0:ppn=14,mem=22gb">> Fusion_eric.sh;
echo "#PBS -l nodes=0:ppn=14,mem=22gb">> Fuscatch_run.sh;
echo "#$tm">> Fusion_eric.sh;
echo "#$tm">> Fuscatch_run.sh;





for i in "${mitico[@]}"
do  pet=(`echo $i|sed 's/_R1_/_R2_/g'`);
pref=(`echo $i|awk -F "/" '{print $NF}'|awk -F "_" '{print $1}' `);
echo "~/software/ericscript-0.4.0/ericscript.pl -name $pref  -p 14 -r /home/sbsuser/databases/bwa_index/hg19.fa -o $nomi/Fusion_res/Eric/"$i" $i $pet;" >> Fusion_eric.sh
done



echo "#!/bin/bash" > Fuscatch_run.sh;
echo "#PBS -l nodes=0:1:ppn=14,mem=22gb">> Fuscatch_run.sh;
echo "#$tm">> Fuscatch_run.sh;

for i in "${mitico[@]}"
do  pet=(`echo $i|sed 's/_R1_/_R2_/g'`);
pref=(`echo $i|awk -F "/" '{print $NF}'|awk -F "_" '{print $1}' `);

  echo "fusioncatcher -d /home/master/database/human_fusioncacther2/ -i $i,$pet  -o $nomi/Fusion_res/Fucatcher/$pref;" >> Fuscatch_run.sh;
    echo "fusioncatcher -d /home/master/database/human_fusioncacther2/ -i $i,$pet  -o  $nomi/Fusion_res/Fucatcher/$pref;" >> Fuscatch_run.sh;

    echo "python2.7  ~/.local/bin/chimerascan_run.py -p 10 --quals solexa ~/database/chimerascan/" $i $pet  $nomi/Fusion_res/chimerascan/$pref >> Chimerascan_run.sh;

done



