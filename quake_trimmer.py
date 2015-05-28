#!/usr/bin/env python

import subprocess
import tempfile
import shutil
import os

def trim(assembler=lambda d: None, reads1=None, reads2=None, trimmer_args="", callback=lambda s, t: None, **kwargs):
	"Not a trimmer per-se; applies QUAKE k-mer error correction to reads. Note that this occurs before, but won't supercede, any correction done internally by the assembler itself."
	if trimmer_args:
		trimmer_args = " " + trimmer_args #pad with a space
	
	temp_path = tempfile.mkdtemp()
	temp_reads1 = os.path.join(temp_path, os.path.basename(reads1).replace(".fastq", ".basic_trimming.fastq"))
	
	
	callback("trimming ({})".format(reads1), status='Running trimming')
	#print "trimming ", temp_reads1, temp_reads2
	subprocess.check_call("quake {}{} {}".format(trimmer_args, reads1, temp_reads1), shell=True)
	kwargs['reads1'] = temp_reads1
	
	if reads2:
		print reads2
		callback("trimming ({})".format(reads2), status='Running trimming')
		temp_reads2 = os.path.join(temp_path, os.path.basename(reads2).replace(".fastq", ".basic_trimming.fastq"))
		subprocess.check_call("quake {}{} {}".format(trimmer_args, reads2, temp_reads2), shell=True)
		kwargs['reads2'] = temp_reads2
	
	
	callback("trimming done ({})".format(temp_path))
	
	
	
	kwargs['callback'] = callback
	
	try:
		d = assembler.assemble(**kwargs)
	
	finally:
		try:
			callback("Cleaning up trimmed reads: {}".format(temp_path))
			shutil.rmtree(temp_path)
		except Exception:
			pass
	
	return d
	
		
if __name__ == "__main__":
	#debug
	import datetime
	import worst_assembler
	def cb(d):
		print "[{}] {}".format(datetime.datetime.today().ctime(), d)
	def bcb(d):
		for (k, v) in d.items():
			cb("{} : {}".format(k, v))
	bcb(trim(path='/home/justin.payne', 
			 assembler=worst_assembler,
			 reads1='/shared/gn2/CFSANgenomes/CFSAN003966/CFSAN003966_01/CFSAN003966_S6_L001_R1_001.fastq', 
			 reads2='/shared/gn2/CFSANgenomes/CFSAN003966/CFSAN003966_01/CFSAN003966_S6_L001_R2_001.fastq', 
			 accession='CFSAN003966_01',
			 insert_size=500,
			 callback=cb,
			 update_callback=bcb,
			 debug=True))