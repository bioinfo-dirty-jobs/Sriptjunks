#!/usr/bin/python
#This is a program  for count Gene on Ensemble data  "process subsection",
#Author: Maurizio Polano, mauriziopolano@blu.it
#Last revision: 03/06/2014

import subprocess
import os,os.path
import sys
import time
import re
from optparse import OptionParser


def main():
	# select inupt
	parser = OptionParser() 
	parser = OptionParser(usage="usage: %prog -i [files  ID,ID,ID]  -d [directory ]  -l list file id -h help  ",
	                          version="%prog 1.0")
	parser.add_option("-d", "--dir", dest="PosDir",type="string",help=" write input file: %prg -i: Please insert the name the Path of the results files [REQUIRED]")
	parser.add_option("-i", "--ID", dest="id_name",type="string",help=" write input file: %prg -i: Please insert the ID files [REQUIRED]")
	parser.add_option("-l", "--list", dest="file_id",type="string",help=" write input file: %prg -i: Please insert the file contain ID")


	(options, args) = parser.parse_args(args=None, values=None)
	if not os.path.exists("Count"):
		os.makedirs("Count")
	if not os.path.exists("Count"):
		os.makedirs("Count")
	if not os.path.exists("de"):
		os.makedirs("de")
    
    
    



	options, args = parser.parse_args()
	
	if options.id_name:

		ID =  options.id_name.split(",")
		
		fileh = open("_count.ense.sh","w")

		for i in ID:
			ctr = options.PosDir+i+"_ens_sort.bam"
			cmd01 =  'samtools sort -no %s de/%s.subset.tmp |samtools view -|/illumina/software/PY276/bin/htseq-count  --mode=intersection-nonempty --stranded=yes --type=exon --idattr=gene_id - /home/sbsuser/databases/bowtie2_ens/Homo_sapiens.GRCh37.72.gtf > Count/%s.count.out'%(ctr,i,i)
			
			print  >> fileh, cmd01

		fileh.close()

	if options.file_id:
		fileh = open("_count.ense.sh","w")
		with open(options.file_id) as p:
			for i in p:
				lines = i.strip("\n")
				ctr = options.PosDir+lines+"_ens_sort.bam"
				cmd01 =  'samtools sort -no %s de/%s.subset.tmp |samtools view -|/illumina/software/PY276/bin/htseq-count  --mode=intersection-nonempty --stranded=yes --type=exon --idattr=gene_id - /home/sbsuser/databases/bowtie2_ens/Homo_sapiens.GRCh37.72.gtf > Count/%s.count.out'%(ctr,lines,lines)
		
				print  >> fileh, cmd01
	fileh.close()

if __name__ == "__main__":
	main()

	
