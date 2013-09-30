#!/usr/bin/env python
import subprocess
import tempfile
import shutil
import os
import re
import sys

import fasta_statter

"""
	
	IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth.

"""
description = "IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth."

core_load = 8 #number of cores this assembler will max out

supports = ('MiSeq', 'IonTorrent')

def assemble(accession, path, reads1, reads2=None, callback=lambda s: None, update_callback=lambda d: None, fasta_file_name=None, debug=True, **kwargs):
	d = {'assembly_version':'',
		 'average_coverage':'',
		 'num_contigs':'',
		 'n50':'',
		 'num_bases':'',
		 'fasta_file':'',
		 'lib_insert_length':'Not determined',
		 'matched':''
		 }
	temp_dir = tempfile.mkdtemp()
	temp_output = tempfile.mkdtemp()
	if not fasta_file_name:
		d['fasta_file'] = fasta_file_name = "{}.idba-ud.fasta".format(accession)
	try:
		#get version
		d['assembly_version'] = 'IDBA-UD v. 1.1.1'
		
		#interleave reads and strip qual scores
		from Bio import SeqIO
		if reads2:
			with open(reads1, 'r') as r1, open(reads2, 'r') as r2, open(os.path.join(temp_dir, 'reads.fasta'), 'w') as reads:
				for (c1, c2) in zip(SeqIO.parse(r1, 'fastq'), SeqIO.parse(r2, 'fastq')):
					reads.write(c1.format('fasta'))
					reads.write(c2.format('fasta'))
		else:
			with open(reads1, 'r') as r1, open(os.path.join(temp_dir, 'reads.fasta'), 'w') as reads:
				for c1 in SeqIO.parse(r1, 'fastq'):
					reads.write(c1.format('fasta'))
		 
		#assemble
		if debug:
			buffer = sys.stdout
		else:
			buffer = open(os.devnull, 'w')
		
		buffer.write(subprocess.check_output("idba -r {} -o {} --pre_correction --num_threads 8".format(os.path.join(temp_dir, 'reads.fasta'), temp_output), shell=True))
		shutil.copyfile(os.path.join(temp_output, 'scaffold.fa'), os.path.join(path, fasta_file_name))
		d.update(fasta_statter.stat_fasta(os.path.join(path, fasta_file_name)))
	
	except subprocess.CalledProcessError as e:
		print type(e), e, e.output
		raise e
	
	finally:
		callback("Cleaning up {}, {}...".format(temp_dir, temp_output))
		#shutil.rmtree(temp_dir)
		#shutil.rmtree(temp_output)
		
	return d
	
	
if __name__ == "__main__":
	#debug
	import datetime
	def cb(d):
		print "[{}] {}".format(datetime.datetime.today().ctime(), d)
	def bcb(d):
		for (k, v) in d.items():
			cb("{} : {}".format(k, v))
	bcb(assemble(path='/home/justin.payne', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R1_001.fastq', 
			 reads2='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R2_001.fastq', 
			 accession='CFSAN001659_01', 
			 callback=cb,
			 update_callback=bcb,
			 debug=True))
	bcb(assemble(path='/home/justin.payne', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R1_001.fastq', 
			 #reads2='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R2_001.fastq', 
			 accession='CFSAN001659_01', 
			 data_type='IonTorrent',
			 callback=cb,
			 update_callback=bcb,
			 debug=True))