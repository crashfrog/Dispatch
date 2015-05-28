#!/usr/bin/env python

import shutil
import os
import glob
import textwrap

from Bio import SeqIO

import fasta_statter

"""
	
	Generic assembler template for writing assembly bindings. Assembly Dispatch provides
	a number of relevant arguments (mostly geared around the requirements of Velvet, but
	this has proved to be robust for most complex de Bruijin assemblers on MiSeq data)
	including the name of a function ("statter") in the fasta_statter module (invoke using
	
		getattr(fasta_statter, statter)(<fasta_file>)
		
	if desired) and a pair of functions for informatory reporting ("callback(string)") and 
	optionally updating the job's information in the database ("update_callback(dict)") by
	dictionary key/value pairs.
	
	Reducing the value of core_load will cause the Dispatch to invoke more jobs, if 
	possible, to maximize CPU usage. The assembler will re-order jobs, if necessary, to
	avoid dispatching more than 8 cores' worth of jobs at a time.
	
	You may assume that your reads have already been trimmed according to the specified
	trimming parameters in the job database; no need to implement trimming here unless
	its assembler-specific. If so, you should call update_callback to reflect this
	trimming.
	
	

"""
description = "WORST: Win On Regular Sequencing Tallies. Cheating assembler. Appends reads to form one big contig until it hits 4.5mb. Worse-case scenario for traditional N50/contig number assembly assessment tools. DO NOT USE"

core_load = 1 #number of cores this assembler will max out

supports = ('MiSeq','IonTorrent')

def assemble(reads1, path, callback=lambda s: None, update_callback=lambda d: None, **kwargs):
	d = {'assembly_version':'NEGATIVE ASSEMBLY CONTROL',
		 'average_coverage':'MISSASSEMBLY DO NOT USE',
		 'num_contigs':'1',
		 'n50':'',
		 'num_bases':'',
		 'fasta_file':'purposeful_misassembly.fasta',
		 'lib_insert_length':'MISASSEMBLY DO NOT USE'
		 }
	try:
		size = 0
		callback("Starting WORST: The Assessment-Testing Purposeful Missassembly.")
		try:
			with open(glob.glob("{}/*.fasta".format(path))[0], 'r') as another_fasta:
				 for c in SeqIO.parse(another_fasta, "fasta"):
					size += len(c)
		except IndexError:
			size = 4500000
		callback("Target assembly size: {} bases.".format(size))
		#assemble
		a = ''
		reads_grabbed = 0
		with open(reads1, 'r') as reads_file:
			reads = SeqIO.parse(reads_file, "fastq")
			while len(a) < size:
				try:
					a += reads.next().seq
				except StopIteration:
					break
				reads_grabbed += 1
				if reads_grabbed % 5000 == 0:
					callback("Appended {} reads.".format(reads_grabbed))
		d['num_bases'] = len(a)
		with open('{}/purposeful_misassembly.fasta'.format(path), 'w') as assembly:
			assembly.write(">PURPOSEFUL_MISASSEMBLY_DO_NOT_USE\n")
			for i in range(0, len(a), 80):
				assembly.write(str(a[i:i+80]))
				assembly.write("\n")
				
		d.update(fasta_statter.stat_fasta('{}/purposeful_misassembly.fasta'.format(path)))
	
	finally:
		pass
		
	return d
	
	
if __name__ == "__main__":
	#debug
	import datetime
	def cb(d):
		print "[{}] {}".format(datetime.datetime.today().ctime(), d)
	def bcb(d):
		for (k, v) in d.items():
			cb("{} : {}".format(k, v))
	print assemble(path='/home/justin.payne', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R1_001.fastq', 
			 reads2='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R2_001.fastq', 
			 accession='CFSAN001659_01', 
			 callback=cb,
			 update_callback=bcb,
			 k_value=144)