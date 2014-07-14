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
import pysam


def main():
	# select inupt
	parser = OptionParser() 
	parser = OptionParser(usage="usage: %prog -i  [file bam ]  -t [gene1,gene2 ] -d [bed file name genes] -e Exons -o [output file] -r hg19 -h help  ",
	                          version="%prog 1.0")
	parser.add_option("-i", "--bam", dest="file_bam",type="string",help=" write input file: %prg -i: Please insert the  bam file [REQUIRED]")
	parser.add_option("-t", "--tag", dest="tag_gene",type="string",help=" write input file: %prg -t: insert gene5' and gene3' [REQUIRED]")
	parser.add_option("-d", "--out", dest="database",type="string",help=" write input file: %prg -i: Please insert the  ouf file")

	parser.add_option("-o", "--outfile", dest="outfile",type="string",help=" write input file: %prg -i: Please insert the  ouf file")
	parser.add_option("-e", "--extr", dest="exon_file",type="string",help=" write input file: %prg -e exon file")
	parser.add_option("-r", "--ref", dest="ref_file",type="string",help=" write input file: %prg -r refernce hg19")
	(options, args) = parser.parse_args(args=None, values=None)
	
    
    



	options, args = parser.parse_args()
	pos = os.path.dirname(os.path.realpath(options.file_bam))

	if options.file_bam:
		print options.file_bam
		gene1,gene2 =  options.tag_gene.split(",")
		gen = {}
		with open(options.database) as p:
			for i in p:
				lines = i.rstrip("\n").split("\t")
				gen.setdefault(lines[3],[]).append(str(lines[0]+":"+lines[1]+"-"+lines[2]))

		if gen.has_key(gene1):
			cord5 = gen[gene1][0]
		if gen.has_key(gene2):
			cord3 = gen[gene2][0]
		print "Start pipeline of finding enconpassing reads on %s -%s "%(gene1,gene2)
		cmd1 = 'samtools view %s %s > %s.txt' % (options.file_bam,cord5,gene1) 
		cmd2 = 'samtools view %s %s > %s.txt' % (options.file_bam,cord3,gene2) 
		fin1 = '%s.txt'% gene1
		fin2 = '%s.txt'% gene2
		print fin2,fin1
		p1 = subprocess.Popen(cmd1,shell=True,stdout=subprocess.PIPE)
		p2 = subprocess.Popen(cmd2,shell=True,stdout=subprocess.PIPE)
		while True:
			if p1.poll() is None:
				time.sleep(3)
				pass
			if p1.poll()==0:
				print 'sam prepared.'
				break
			if p1.poll() is not None and p1.poll() != 0:
				raise Exception('Error building junction db index')

		while True:
			if p2.poll() is None:
				time.sleep(3)
				pass
			if p2.poll()==0:
				print 'sam prepared.'
				break
			if p2.poll() is not None and p2.poll() != 0:
				raise Exception('Error building samtools')
		t = 1
		print "Write out files of encompassing reads and split reads %s -%s "%(gene1,gene2)
		fileh = open(options.outfile,"w")
		fileh2 = open(pos+"/"+"reads.name.map","w")
		with open(pos+'/'+fin1) as p:
			for i in p:
				lines = i.strip("\n").split("\t")
				
				if lines[2] != lines[6]:
					if not (lines[6]  == "=" or lines[6] == "*"):
						ct2 = ">"+lines[0]+"_"+lines[2]+lines[6]+"\n"+lines[9]
						ct = ">"+str(t)+ "_"+lines[2]+lines[6]+"\n"+lines[9]
						print  >> fileh, ct
						print >> fileh2, ct2
				t =t+1
		with open(pos+'/'+fin2) as p:
			for i in p:
				lines = i.strip("\n").split("\t")
				
				if lines[2] != lines[6]:
					if not (lines[6]  == "=" or lines[6] == "*"):
						ct2 = ">"+lines[0]+"_"+lines[2]+lines[6]+"\n"+lines[9]
						
						ct = ">"+str(t)+ "_"+lines[2]+lines[6]+"\n"+lines[9]
						print  >> fileh, ct
						print >> fileh2, ct2
				t =t+1
		fileh.close()
		os.remove(fin1)
		os.remove(fin2)
		fileh2.close()
		

		#prepare exon database
		exon = {}
		pot = os.path.dirname(os.path.realpath(options.exon_file))
		pot2 = os.path.dirname(os.path.realpath(options.ref_file))

		with open(options.exon_file) as p:
			for i in p:
				lines = i.rstrip("\n").split(" ")
				posi = lines[0]+":"+lines[1]+"-"+lines[2]
				gene = lines[3]+"_"+lines[5]
				exon.setdefault(gene,[]).append(posi)

			
		for i in exon.keys():

			if i.startswith(gene1):
				cos=  exon[i][0]
				cmd4 = 'samtools faidx   %s  %s >> de.tmp ' % (options.ref_file,cos)
				p4 = subprocess.Popen(cmd4,shell=True,stdout=subprocess.PIPE)
				while True:
					if p4.poll() is None:
						time.sleep(2)
						pass
					if p4.poll()==0:
						print 'Create fasta file of %s' %(gene1)
						break
					if p4.poll() is not None and p4.poll() != 0:
						raise Exception('Error on faidx')


			if i.startswith(gene2):
				cos=  exon[i][0]
				print exon[i][0],i
				cmd5 = 'samtools faidx   %s  %s >> de.tmp ' % (options.ref_file,cos)
				p5 = subprocess.Popen(cmd5,shell=True,stdout=subprocess.PIPE)
				while True:
					if p5.poll() is None:
						time.sleep(2)
						pass
					if p5.poll()==0:
						print 'Create fasta file of %s' %(gene2)
						break
					if p5.poll() is not None and p5.poll() != 0:
						raise Exception('Error on faidx')



		#CREATE DATABASE FOR BLAT
		print "Create database of exon of  %s -%s "%(gene1,gene2)

		tmp = pos+'/'+'de.tmp'
		cmd3 = 'faToTwoBit -noMask  %s  %s.2bit -ignoreDups' % (tmp,'gene_gene')
		p3 = subprocess.Popen(cmd3,shell=True,stdout=subprocess.PIPE)

		while True:
			if p1.poll() is None:
				time.sleep(3)
				pass
			if p1.poll()==0:
				print '#'
				break
			if p1.poll() is not None and p1.poll() != 0:
				raise Exception('Error building junction gene-gene')

		
		print "Start Blat of split reads..."
		n = pos+"/"+'gene_gene.2bit'
		name4= pos+"/"+"reads.name.map"
		name5 =pos+"/"+options.outfile
		final_results = '%s_%s-%s' %(gene1,gene2,"gene_gene")
		hy ='blat -t=DNA -q=RNA -stepSize=5 -tileSize=11  %s  %s  %s.psl'%(n,name5,final_results)
		print hy
		p5 = subprocess.Popen(hy,shell=True,stdout=subprocess.PIPE)

		while True:
			if p5.poll() is None:
				time.sleep(3)
				pass
			if p5.poll()==0:
				print '#'
				break
			if p5.poll() is not None and p5.poll() != 0:
				raise Exception('Error running blat')

		os.remove("de.tmp")



if __name__ == "__main__":
	main()
		