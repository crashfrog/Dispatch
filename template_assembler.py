import subprocess
import tempfile
import shutil
import os
import re

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
description = "One sentence description of the assembler being invoked."

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