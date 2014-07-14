#!/bin/bash


ID=$1;
mkdir report
mkdir tmp
mkdir OUTPUT
#sort ddati
samtools sort accepted_hits.bam $ID"_sort"; 
#add RG group
java -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/ -jar  -Xmx5g /illumina/software/PROG/picard-tools-1.99/AddOrReplaceReadGroups.jar I=$ID"_sort.bam" O=$ID"_RG.bam" RGPL=illumina RGSM=12358 RGPU=ID1324 LB=1324 ID=$ID;
#Sort_data
java -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/ -jar -Xmx5g /illumina/software/PROG/picard-tools-1.99/SortSam.jar   INPUT=$ID"_RG.bam"  OUTPUT=$ID"_RG_sorted_reads.bam"  SORT_ORDER=coordinate  VALIDATION_STRINGENCY=LENIENT TMP_DIR=/home/sbsuser/mauriziotemp/ ;

#MarkDuplicates
java -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/ -jar -Xmx5g  /illumina/software/PROG/picard-tools-1.99/MarkDuplicates.jar I=$ID"_RG_sorted_reads.bam"  O=$ID"_RG_sorted_reads_dp.bam" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics;
#Sort_data
java -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/ -jar -Xmx5g /illumina/software/PROG/picard-tools-1.99/SortSam.jar   INPUT=$ID"_RG_sorted_reads_dp.bam"  OUTPUT=$ID"_RG_sorted_reads_dp_sort2.bam"  SORT_ORDER=coordinate  VALIDATION_STRINGENCY=LENIENT TMP_DIR=/home/sbsuser/mauriziotemp/;

#Reorder using ucsf
java -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/ -jar -Xmx5g  -jar /illumina/software/PROG/picard-tools-1.99/ReorderSam.jar INPUT=$ID"_RG_sorted_reads_dp_sort2.bam"  OUTPUT=$ID"_RG_sorted_reads_dp_sort2_ucs.bam" REFERENCE=/illumina/software/database/ucsc.hg19.fasta



#index the file
java -Djava.io.tmpdir=tmp/ -jar -Xmx5g  -jar /illumina/software/PROG/picard-tools-1.99/BuildBamIndex.jar INPUT=$ID"_RG_sorted_reads_dp_sort2_ucs.bam"
#SplitNCigarRead
java  -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/  -Xmx10g -jar /illumina/software/PROG/GATK31/GenomeAnalysisTK.jar -T SplitNCigarReads -R /illumina/software/database/ucsc.hg19.fasta -I $ID"_RG_sorted_reads_dp_sort2_ucs.bam" -o $ID"_RG_sorted_reads_dp_sort2_ucs_split.bam" -rf ReassignMappingQuality -DMQ 60 -U ALLOW_N_CIGAR_READS

#Reorder sam
java -Djava.io.tmpdir=/illumina/runs/FASTQ/Analisidic2013/ALTRO2/risultati/Mutazioni/138/tmp/ -Xmx10g -jar /illumina/software/PROG/picard-tools-1.99/ReorderSam.jar INPUT=$ID"_RG_sorted_reads_dp_sort2_ucs_split.bam"  OUTPUT=$ID"_RG_sorted_reads_dp_sort2_ucs_split_ord.bam" REFERENCE=/illumina/software/database/ucsc.hg19.fasta

#base recalibration
java  -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/   -Xmx10g -jar /illumina/software/PROG/GATK31/GenomeAnalysisTK.jar -T BaseRecalibrator  -R /illumina/software/database/ucsc.hg19.fasta  --knownSites /illumina/software/database/dbsnp_138.hg19.excluding_sites_after_129.vcf -I $ID"_RG_sorted_reads_dp_sort2_ucs_split.bam" -o report/$ID"__report.grp" -l INFO 


#base build index
java -Djava.io.tmpdir=tmp/ -jar -Xmx5g  -jar /illumina/software/PROG/picard-tools-1.99/BuildBamIndex.jar INPUT=$ID"_RG_sorted_reads_dp_sort2_ucs_split_ord.bam"


#togli le reads non buone

java  -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/  -Xmx10g -jar /illumina/software/PROG/GATK31/GenomeAnalysisTK.jar -T PrintReads -R /illumina/software/database/ucsc.hg19.fasta -I $ID"_RG_sorted_reads_dp_sort2_ucs_split_ord.bam" -BQSR report/$ID"__report.grp" -o $ID"_RG_sorted_reads_dp_sort2_ucs_split_ord_recal.bam"


# fai la chiamata
java  -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/  -Xmx5g -jar /illumina/software/PROG/GATK31/GenomeAnalysisTK.jar -T HaplotypeCaller -R /illumina/software/database/ucsc.hg19.fasta -I $ID"_RG_sorted_reads_dp_sort2_ucs_split_ord_recal.bam"  -recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o out/$ID".vcf" -L  /home/sbsuser/databases/Target_gene/GIST.target.bed


#Filtra via
java  -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/  -Xmx5g -jar /illumina/software/PROG/GATK31/GenomeAnalysisTK.jar -T VariantFiltration -R /illumina/software/database/ucsc.hg19.fasta -V  out/$ID".vcf"  -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o out/$ID"_filtrato.vcf"

#converti in vcf4
/illumina/software/PROG/annovar2/convert2annovar.pl  --format vcf4 out/$ID"_filtrato.vcf" > out/$ID"_filtrato.avinput";


qsub -cwd -v PATH -b y /illumina/software/PROG/annovar2/table_annovar.pl out/$ID"_filtrato.avinput"  /illumina/software/PROG/annovar/humandb/ -buildver hg19 -out out/$ID"_mianno.xls" -remove -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp138,ljb2_all -operation g,r,r,f,f,f,f -nastring NA -csvout 
