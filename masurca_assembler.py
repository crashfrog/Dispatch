#!/usr/bin/env python

import subprocess
import tempfile
import shutil
import os

"""

	MaSuRCA assembler binding for Assembly Dispatch.
	

"""

manifests = {'MiSeq':"""
DATA
PE= pe 241 20 {reads1} {reads2}
END

PARAMETERS
#this is k-mer size for deBruijn graph values between 25 and 101 are supported, auto will compute the optimal size based on the read data and GC content. Do not set this longer than PE read length!!!
GRAPH_KMER_SIZE=auto
#set this to 1 for Illumina-only assemblies and to 0 if you have 2x or more long (Sanger, 454) reads
USE_LINKING_MATES=1
#this parameter is useful if you have too many jumping library mates. Typically set it to 60 for bacteria and something large (300) for mammals
#LIMIT_JUMP_COVERAGE = 60
#these are the additional parameters to Celera Assembler. do not worry about performance, number or processors or batch sizes -- these are computed automatically. for mammals do not set cgwErrorRate above 0.15!!!
CA_PARAMETERS = ovlMerSize=30 cgwErrorRate=0.25 ovlMemory=4GB
#auto-detected number of cpus to use
NUM_THREADS=8
#this is mandatory jellyfish hash size 10x the genome size is a good starting value
JF_SIZE=45000000
#this specifies if we do (1) or do not (0) want to trim long runs of homopolymers (e.g. GGGGGGGG) from 3' read ends, use it for high GC genomes
DO_HOMOPOLYMER_TRIM=0
END	
"""}

description = "MaSuRCA, an experimental combination OLC/de Brujin assembler developed at UMD. Not yet supported."

core_load = 8

supports = manifests.keys()

def assemble(data_type='MiSeq', callback=lambda s: None, update_callback=lambda d: None, **kwargs):
	
	#callback("MaSuRCA not currently supported.")
	#raise ValueError("MaSuRCA not currently supported.")
	
	temp_path = tempfile.mkdtemp()
	
#write MaSuRCA config file
	with open(os.path.join(temp_path, "{accession}_config.txt".format(**kwargs)), 'w') as config_file:
		config_file.write("""PATHS
JELLYFISH_PATH=/install_path/bin
SR_PATH=/install_path/bin
CA_PATH=/install_path/CA/Linux-amd64/bin
END

#DATA is specified as type {PE,JUMP,OTHER}= two_letter_prefix mean stdev fastq(.gz)_fwd_reads fastq(.gz)_rev_reads
#NOTE that PE reads are always assumed to be innies, i.e. ---> <---, and JUMP are assumed to be outties <--- --->; if there are any jump libraries that are innies, such as longjump, specify them as JUMP and specify NEGATIVE mean
#IT IS MANDATORY to supply some paired end reads
#reverse (R2) reads are optional for PE libraries and mandatory for JUMP libraries
#any OTHER sequence data (454, Sanger, Ion torrent, etc) must be first converted into Celera Assembler compatible .frg files (see http://wgsassembler.sourceforge.com)

{params}
""".format(params=manifests[data_type].format(**kwargs), **kwargs))

	try:
		subprocess.check_call(
	finally:
		print "Cleaning up {}...".format(temp_path)
		shutil.rmtree(temp_path)





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
			 k_value=144))