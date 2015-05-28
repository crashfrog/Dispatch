#!/usr/bin/env python

import subprocess
import tempfile
import shutil
import os
import re


import fasta_statter

"""
	
	Assembly Dispatch bindings for SPAdes Single-Cell assembler (Bankevitch et al,
	J Comput Biol, 2012), a de Brujin graph assembler for paired-end single-cell
	sequencing.
	
	Trying something new here - using file-like object to try to parse SPAdes output as it 
	runs. Also SPAdes has a Python script wrapper, so this one is a little different - we 
	actually import the SPAdes module and invoke its main method with an arg vector.
	
	Jan 6 2014: No, we don't. That was stupid from the get-go.
	

"""

class SpadesParser(object):
	"File-like object to parse IO streams and invoke a function in response to a regex match."

	def __init__(self, function, regex, buffer):
		self.func = function
		self.re = re.compile(regex, re.MULTILINE)
		self.buffer = buffer
		self.passthrough = False
		
	def write(self, s):
		"Capture a text stream, look for regex matches; if found call func with the match"
		if 'Error' in s:
			self.passthrough=True
		r = self.re.search(s)
		if r and not self.passthrough:
			self.func(r.groups()[0])
		if self.passthrough:
			self.func(s)
		return self.buffer.write(s)
		
	def fileno(self):
		return 0
		
	def flush(self):
		pass
		
	def close(self):
		pass
		
		
description = "SPAdes, a de Brujin graph assembler with advanced error correction and handling of paired-end data."

core_load = 1 #number of cores this assembler will max out

#spades_exec = '/home/justin.payne/SPAdes-2.5.0-Linux/bin/'
#spades_exec = '/home/justin.payne/SPAdes-3.0.0-Linux/bin/'
spades_exec = ''

supports = ('MiSeq', 'IonTorrent')

def assemble(path, accession, fasta_file_name=None, callback=lambda s: None, update_callback=lambda d: None, debug=True, **kwargs):
	
	import sys
	if not fasta_file_name:
		fasta_file_name = "{}.spades.fasta".format(accession)
	kwargs['temp_dir'] = tempfile.mkdtemp()
	
	d = {'assembly_version':'SPAdes v. 3.0.0',
		 'average_coverage':'',
		 'num_contigs':'',
		 'n50':'',
		 'num_bases':'',
		 'fasta_file':fasta_file_name,
		 'lib_insert_length':'Not determined',
		 'matched':'-'
		 }
		 

	try:
		#import the spades python wrapper
		#sys.path.append(spades_exec)
		#import spades
		#get version
		#d['assembly_version'] = "SPAdes v. {}".format(spades.spades_version.replace("\n", ""))
		
		#assemble
		#print "version:", d['assembly_version']
		callback("running SPAdes in {temp_dir}...".format(**kwargs))
		
		def status_callback(r):
			callback("running: " + r)
		
		if not debug:
			sys.stdout = SpadesParser(status_callback, r"^ *=+ ([A-Za-z0-9;:., ]+)", open(os.devnull, 'w'))
		#argv = "--disable-gzip-output --careful -t 8 -m 64 -1 {reads1} -2 {reads2} -o {temp_dir}".format(**kwargs).split(" ")
		#if debug:
		#	callback("spades {}".format(" ".join(argv)))
		#spades.main(argv)
		
		if 'reads2' in kwargs and kwargs['reads2']:
			subprocess.check_call("{spades_exec}spades -t 8 -m 64 -1 {reads1} -2 {reads2} -o {temp_dir}".format(spades_exec=spades_exec, **kwargs), shell=True)
		else:
			subprocess.check_call("{spades_exec}spades --iontorrent -t 8 -m 64 -s {reads1} -o {temp_dir}".format(spades_exec=spades_exec, **kwargs), shell=True)
		
		update_callback({"k_value":'21, 33, 55, 77, 99, 127'})
		
		callback("Copying {temp_dir}/contigs.fasta...".format(**kwargs))
		shutil.copyfile("{temp_dir}/contigs.fasta".format(**kwargs), "{}/{}".format(path, fasta_file_name))
		shutil.copyfile("{temp_dir}/spades.log".format(**kwargs), "{}/spades.log".format(path))

		d.update(fasta_statter.stat_velvet(os.path.join(path, fasta_file_name), 33))
		
#	except subprocess.CalledProcessError as e:
#		raise ValueError(str(type(e)) + str(e) + str(e.output))
		
	except Exception as e:
		if debug:
			import traceback
			traceback.print_exc(sys.stdout)
			raise e
		raise ValueError("SPAdes assembly failure.")
			
	
	finally:
		callback("Cleaning up {temp_dir}...".format(**kwargs))
		shutil.rmtree(kwargs['temp_dir'])
		sys.stdout = sys.__stdout__
		
	return d
	
	
if __name__ == "__main__":
	import sys
	import datetime
	def cb(s):
		print "[{}]".format(datetime.datetime.today().ctime()), s
	def ucb(d):
		for (key, value) in d.items():
			print key, ":", value
	print assemble(path='/home/justin.payne/', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001656/CFSAN001656_01/CFSAN001656_S8_L001_R1_001.fastq', 
			 reads2='/shared/gn2/CFSANgenomes/CFSAN001656/CFSAN001656_01/CFSAN001656_S8_L001_R2_001.fastq', 
			 accession='CFSAN001656_01', 
			 callback=cb,
			 update_callback=ucb,
			 debug=True)
	print assemble(path='/home/justin.payne', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN006329/CFSAN006329_01/CFSAN006329.reads.fastq', 
			 accession='CFSAN006329_01', 
			 callback=cb,
			 debug=True,
			 update_callback=ucb)