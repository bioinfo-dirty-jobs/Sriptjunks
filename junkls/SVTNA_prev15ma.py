#!/usr/bin/env python
#This is a program to implement the RNA_calling  pipeline "process subsection",
#Author: Maurizio Polano, mauriziopolano@blu.it
#Last revision: 22/04/2014

import subprocess
import os,os.path
import sys
import time
import re
import random
args=sys.argv

help_menu='''\nPipeline for Mutation call on RNAseq Data SVRNA.
\n\t**Usage**:
\tcode -sample 
\t**Parameters**:
\t-sample		name complete of sample to use
\t -inputdir if the data are not in the same position of script
\t- ID		ID identify the analisys step 
\t -step start with 01
\t-h		print help message
'''


if not os.path.exists("report"):
    os.makedirs("report")
if not os.path.exists("out"):
    os.makedirs("out")


if not os.path.exists("Finale"):
    os.makedirs("Finale")

param = "'$4<6 '"
if '-h' in args or '-help' in args or len(args)==1:
    print help_menu
    sys.exit(0)
########################################################################################
#underlying utilities, automatically detected

#Default uses 12 nodes in HPC
#########################################################################################
if '-sample' not in args:
    sys.exit('ERROR: Sample name is needed')

i=args.index('-sample')
sample=args[i+1]
if '-inputdir' not in args:
    inputpath=os.path.abspath('./')
else:
    i=args.index('-inputdir')
    inputpath=args[i+1]


if '-step_info' in args:
    print steps_info
    sys.exit(0)


bampath=os.path.abspath(inputpath)+'/%s'%sample
i=args.index('-step')
step=args[i+1]


####################################################################################################


LB= random.randrange(1,2000)
i=args.index('-ID')
if '-ID' in args:
    ID = str(args[i+1])
else:
    ID= sample


RGPU=sample+str(LB)

########################################################################################
#pipeline command lines.

step_01_cmd=['samtools sort  -m 1000000000 %s %s.sorted.bam'%(sample,ID)]
step_02_cmd=['java -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/ -jar   /illumina/software/PROG/picard-tools-1.99/AddOrReplaceReadGroups.jar I=%s.sorted.bam.bam O=%s_RG.bam RGPL=illumina RGSM=%s RGPU=%s LB=%s ID=%s;'%(ID,ID,LB,RGPU,LB,ID)]
step_03_cmd=['samtools fixmate -r %s_RG.bam %s_RG_fixmate.bam'%(ID,ID)]
step_04_cmd=['java  -Xmx5g  -XX:+UseConcMarkSweepGC -XX:-UseGCOverheadLimit  -jar  /illumina/software/PROG/picard-tools-1.99/SortSam.jar   INPUT=%s_RG_fixmate.bam  OUTPUT=%s_RG_sorted_reads.bam  SORT_ORDER=coordinate  VALIDATION_STRINGENCY=LENIENT TMP_DIR=/home/sbsuser/maurizio2tmp/ MAX_RECORDS_IN_RAM=1000000'%(ID,ID)]
step_05_cmd=['java -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/ -Xmx20g -jar   /illumina/software/PROG/picard-tools-1.99/MarkDuplicates.jar I=%s_RG_sorted_reads.bam  O=%s_RG_sorted_reads_dp.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=true'%(ID,ID)]
step_06_cmd=['java -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/ -jar  /illumina/software/PROG/picard-tools-1.99/SortSam.jar   INPUT=%s_RG_sorted_reads_dp.bam  OUTPUT=%s_RG_sorted_reads_dp_sort2.bam  SORT_ORDER=coordinate  VALIDATION_STRINGENCY=LENIENT TMP_DIR=/home/sbsuser/mauriziotemp/'%(ID,ID)]

step_07_cmd=['java -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/ -jar -Xmx20g -jar /illumina/software/PROG/picard-tools-1.99/ReorderSam.jar INPUT=%s_RG_sorted_reads_dp_sort2.bam  OUTPUT=%s_RG_sorted_reads_dp_sort2_ucs.bam REFERENCE=/illumina/software/database/ucsc.hg19.fasta'%(ID,ID)]
step_08_cmd=['java -Djava.io.tmpdir=tmp/ -jar -Xmx5g  -jar /illumina/software/PROG/picard-tools-1.99/BuildBamIndex.jar INPUT=%s_RG_sorted_reads_dp_sort2_ucs.bam'%(ID)]
step_09_cmd=['java  -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/  -Xmx10g -jar /illumina/software/PROG/GATK31/GenomeAnalysisTK.jar -T SplitNCigarReads -R /illumina/software/database/ucsc.hg19.fasta -I %s_RG_sorted_reads_dp_sort2_ucs.bam -o %s_RG_sorted_reads_dp_sort2_ucs_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS'%(ID,ID)]
step_10_cmd=['java -Djava.io.tmpdir=/illumina/runs/FASTQ/Analisidic2013/ALTRO2/risultati/Mutazioni/138/tmp/ -Xmx10g -jar /illumina/software/PROG/picard-tools-1.99/ReorderSam.jar INPUT=%s_RG_sorted_reads_dp_sort2_ucs_split.bam  OUTPUT=%s_RG_sorted_reads_dp_sort2_ucs_split_ord.bam REFERENCE=/illumina/software/database/ucsc.hg19.fasta VALIDATION_STRINGENCY=SILENT '%(ID,ID)]
step_11_cmd=['java  -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/   -Xmx10g -jar /illumina/software/PROG/GATK31/GenomeAnalysisTK.jar -T BaseRecalibrator  -R /illumina/software/database/ucsc.hg19.fasta  --knownSites /illumina/software/database/dbsnp_138.hg19.excluding_sites_after_129.vcf -I %s_RG_sorted_reads_dp_sort2_ucs_split.bam -o report/%s__report.grp -l INFO '%(ID,ID)]
step_12_cmd=['java -Djava.io.tmpdir=tmp/ -jar -Xmx5g  -jar /illumina/software/PROG/picard-tools-1.99/BuildBamIndex.jar INPUT=%s_RG_sorted_reads_dp_sort2_ucs_split_ord.bam'%(ID)]
step_13_cmd=['java  -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/  -Xmx10g -jar /illumina/software/PROG/GATK31/GenomeAnalysisTK.jar -T PrintReads -R /illumina/software/database/ucsc.hg19.fasta -I %s_RG_sorted_reads_dp_sort2_ucs_split_ord.bam -BQSR report/%s__report.grp -o %s_RG_sorted_reads_dp_sort2_ucs_split_ord_recal.bam'%(ID,ID,ID)]
step_14_cmd=['java  -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/  -Xmx5g -jar /illumina/software/PROG/GATK31/GenomeAnalysisTK.jar -T HaplotypeCaller -R /illumina/software/database/ucsc.hg19.fasta -I %s_RG_sorted_reads_dp_sort2_ucs_split_ord_recal.bam  -recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o out/%s.vcf -L  /home/sbsuser/databases/Target_gene/GIST.target.bed'%(ID,ID)]
step_15_cmd=['java  -Djava.io.tmpdir=/home/sbsuser/mauriziotemp/  -Xmx5g -jar /illumina/software/PROG/GATK31/GenomeAnalysisTK.jar -T VariantFiltration -R /illumina/software/database/ucsc.hg19.fasta -V  out/%s.vcf  -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o out/%s_filtrato.vcf'%(ID,ID)]
step_16_cmd=['/illumina/software/PROG/annovar2/convert2annovar.pl  --format vcf4 out/%s_filtrato.vcf > out/%s_filtrato.avinput'%(ID,ID)]
step_17_cmd=['qsub -cwd -v PATH -b y /illumina/software/PROG/annovar2/table_annovar.pl out/%s_filtrato.avinput  /illumina/software/PROG/annovar/humandb/ -buildver hg19 -out out/%s_mianno.xls -remove -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp138,ljb2_all -operation g,r,r,f,f,f,f -nastring NA -csvout'%(ID,ID)]
step_18_cmd=['bedtools genomecov -ibam %s_RG_sorted_reads_dp_sort2_ucs_split_ord_recal.bam  -bga -split|awk %s |intersectBed  -a ~/databases/Target_gene/GIST.target.bed -b - > out/%s_zone_lowcoverage.bed'%(ID,param,ID)]
#post_6_2_clean=['rm -f %s.withRG.paired.sorted.bam'%sample,'rm -f %s.withRG.paired.sorted.bam.bai'%sample,'rm -f %s.orig.csv'%sample]
step_19_cmd=['rm  %s_RG_sorted_reads_dp_sort2_ucs_split.bam %s_RG_sorted_reads_dp_sort2_ucs.bam %s_RG.bam  %s_RG_sorted_reads_dp_sort2_ucs_split_ord_recal.bai'%(ID,ID,ID,ID)]
step_20_cmd=['/illumina/software/PY276/bin/bam_stat.py -i %s.sorted.bam.bam   &> out/%sstatistics.txt'%(sample,ID)]

step_21_cmd=['/illumina/software/PY276/bin/junction_annotation.py -i %s.sorted.bam.bam -r /illumina/software/database/RSeqc/hg19_RefSeq.bed  -o %s &> out/%s_annotaion.juction.txt'%(sample,ID,ID)]

step_22_cmd=['mv out/%s_annotaion.juction.txt Finale/'%(ID)]
step_23_cmd=['mv out/%sstatistics.txt  Finale/'%(ID)]
step_24_cmd=['mv out/%s Finale/'%(ID)]
step_25_cmd=['samtools index %s_RG_sorted_reads_dp_sort2_ucs_split_ord_recal.bam '%(ID)]
step_26_cmd=['intersectBed -abam %s_RG_sorted_reads_dp_sort2_ucs_split_ord_recal.bam -b ~/databases/Target_gene/GIST.target.bed > Finale/%s_Gist_target.bam'%(ID,ID)]
step_27_cmd=['samtools index Finale/%s_Gist_target.bam'%(ID)]
step_28_cmd=['tar cvfz %s_tar.gz  Finale/'%(ID)]



#step_7_cmd=['java -Xmx8g -jar %s/MarkDuplicates.jar I=%s.withRG.GATKRecalibrated.bam O=%s.withRG.GATKRecalibrated.flagged.bam METRICS_FILE=%s.Duplicates_metrics.txt VALIDATION_STRINGENCY=SILENT TMP_DIR=tmp/'%(picard,sample,sample,sample),'%s index %s.withRG.GATKRecalibrated.flagged.bam'%(samtools,sample)]
#post_7_clean=['rm -f %s.withRG.GATKRecalibrated.bam'%sample,'rm -f %s.withRG.GATKRecalibrated.bai'%sample,'rm -f %s.Duplicates_metrics.txt'%sample]
#step_8_cmd=["java -Xmx8g -jar %s -ttype 2 -t %s -r %s -s '%s|%s.withRG.GATKRecalibrated.flagged.bam|Disc' -o %s/"%(seqc,genome_gtf,genome_fasta,sample,sample,sample)]
#post_8_clean=[]

cmdset=[]
for item in globals().keys():
    if item.endswith('_cmd'):
        cmdset.append(item)
cmdset.sort()

print cmdset



try:
    cmd_entry=cmdset.index('step_'+step+'_cmd')    ##entry point
except ValueError:
    sys.exit('ERROR: STEP not recognized')

def _parsecmd(cmd):
    '''pase cmd str for step information'''
    info=cmd.split('_')
    cstep='_'.join(info[1:-1])
    return cstep

fine =inputpath +"/Script"+sample+ID+".sh"
###########################
##headers
outfile=open(fine,'w')
outfile.write('#! /bin/sh\n')


#commands
outfile.write('echo "Job start: `date`"\n')

for i in range(0,len(cmdset)):
    cmd=cmdset[i]
    
    clean_flag=1
    outfile.write('echo "step %s start: `date`"\n'%i)
    outfile.write('if\n')
    outfile.write('\t%s\n'%('\n'.join(eval(cmd))))
    outfile.write('then\n')
    outfile.write('\techo "step %s done: `date`"\n'%i)
    outfile.write('else\n')
    outfile.write('\techo "step %s ERROR"\n'%i)
    outfile.write('\texit\n')
    outfile.write('fi\n')

   
    

outfile.write('echo "PIPELINE FINISHED"\n')
outfile.close()
######################################################################
