import subprocess
import tempfile
import shutil
import os
import re

import fasta_statter

"""
	
	SGA - String Graph Assembler
	
	Jared T. Simpson and Richard Durbin
	Efficient construction of an assembly string graph using the FM-index.
	Bioinformatics (2010) 26 (12): i367-i373 doi:10.1093/bioinformatics/btq217
	
	

"""
description = "String graph assembler for microbial reads. Uses FM-index and BW-transform to efficiently find overlaps between reads."

core_load = 8 #number of cores this assembler will max out

supports = ('MiSeq','IonTorrent','PacBio','454')

def assemble(callback=lambda s: None, update_callback=lambda d: None, k_value, **kwargs)
	d = {'assembly_version':'',
		 'average_coverage':'',
		 'num_contigs':'',
		 'n50':'',
		 'num_bases':'',
		 'fasta_file':'',
		 'lib_insert_length':'Not determined',
		 'matched':''
		 }
	try:
		#get version
		 
		#assemble
	
	d.update(fasta_statter.stat_fasta(...))
	
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
	print assemble(path='/home/justin.payne', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001812/CFSAN001812_01/R_2012_04_12_16_25_14_user_IT1-33-BW3_Auto_IT1-33-BW3_38.fastq', 
			 accession='CFSAN001812_01', 
			 callback=cb,
			 update_callback=bcb,
			 data_type='IonTorrent')