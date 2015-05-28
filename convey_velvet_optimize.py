#!/usr/bin/env python

import subprocess
import shutil
import os
import tempfile
import datetime
import fasta_statter


import convey_velvet_assembler
import velvet_optimize

"""
	
	Wrapper to execute Velvet_optimize around Convey's implementation of velvet.

"""

description = "Two-round local-maximum k-value optimization assembly using Convey's implementation of Velvet."

core_load = 1

supports = velvet_optimize.supports #supports whatever the underlying Velvet binding supports

def assemble(*args, **kwargs):
	velvet_optimize.velvet_assembler_notrimmer = convey_assembler
	velvet_optimize.assemble(*args, **kwargs)
	
	
if __name__ == "__main__":
	#debug
	import datetime
	def cb(d):
		print "[{}] {}".format(datetime.datetime.today().ctime(), d)
	def bcb(d):
		for (k, v) in d.items():
			cb("{} : {}".format(k, v))
	bcb(assemble(path='/shared/gn2/CFSANgenomes/CFSAN001659/asm/', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R1_001.fastq', 
			 reads2='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R2_001.fastq', 
			 accession='CFSAN001659_01', 
			 callback=cb,
			 update_callback=bcb,
			 k_value=144,
			 debug=True))
	bcb(assemble(path='/shared/gn2/CFSANgenomes/CFSAN001659/asm/', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R1_001.fastq', 
			 #reads2='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R2_001.fastq', 
			 accession='CFSAN001659_01', 
			 data_type='IonTorrent',
			 callback=cb,
			 update_callback=bcb,
			 k_value=144,
			 debug=True))