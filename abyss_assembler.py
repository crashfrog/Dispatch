#!/usr/bin/env python

import subprocess
import shutil
import os
import re
import tempfile
import fasta_statter

"""
	Assembly Dispatch binding for ABySS, the Assembly By Short Sequences assembler
	(Simpson et al Genome Research 2009) (http://www.bcgsc.ca/platform/bioinfo/software/abyss)

"""

core_load = 8 #number of CPU cores this assembler will max out

description = "ABySS, a fast de Brujin assembler for extremely short reads."

supports = ('MiSeq',)

def assemble(path, fasta_file_name=None, callback=lambda s: None, update_callback=lambda d: None, debug=False, **parameters):

	d = {'assembly_version':'',
		 'average_coverage':'',
		 'num_contigs':'',
		 'n50':'',
		 'num_bases':'',
		 'fasta_file':'',
		 'lib_insert_length':'Not determined',
		 'matched':''
		 }
		 
	if not debug:
		parameters['debug'] = " 2> /dev/null"
	else:
		parameters['debug'] = ''
	
	if int(parameters['k_value']) > 64:
		update_callback({'k_value':64})
		 
	parameters['k_value'] = max(min(int(parameters['k_value']), 64), 16)
	

	if not fasta_file_name:
		fasta_file_name = "{}.fasta".format(parameters['accession'])
	d['fasta_file'] = fasta_file_name
		
	parameters['working_dir'] = tempfile.mkdtemp() #get a temp dir for working files
	
	try:	
		#make an abyss output file, starting with the CLI invocation
		results = "abyss-pe np=8 k={k_value} name={accession} in='{reads1} {reads2}' -C {working_dir}\n".format(**parameters) 
		#raise ValueError('ABySS not yet implemented.')
		callback('running abyss-pe in {working_dir}'.format(**parameters))
		results += subprocess.check_output("abyss-pe np=8 k={k_value} name={accession} in='{reads1} {reads2}' -C {working_dir}{debug}".format(**parameters), shell=True)
	
	
		m = re.search(r"ABySS (\d*.\d*.\d*)", results)
		if m:
			d['assembly_version'] = "ABySS v. {}".format(m.group(1))
		
		
		callback("Statting assembly...")
		d.update(fasta_statter.stat_abyss("{working_dir}/{accession}-contigs.fa".format(**parameters), k_value=parameters['k_value']))
		
		shutil.copyfile("{working_dir}/{accession}-contigs.fa".format(**parameters), os.path.join(path, fasta_file_name))
		
		d['fasta_file'] = fasta_file_name
		
		#results += fasta_statter.read_map(os.path.join(path, fasta_file_name), parameters['reads1'], parameters['reads2'])
	
	finally:
		try:
			callback("Cleaning up {working_dir}...".format(**parameters))
			shutil.rmtree(parameters['working_dir'])
			with open(os.path.join(path, "abyss.stats"), 'w') as statsfile:
				statsfile.write(results)
		except OSError:
			pass
		except IOError:
			pass
		
	return d
	
	
if __name__ == "__main__":
	#debug
	import datetime
	import sys
	def cb(d):
		print "[{}] {}".format(datetime.datetime.today().ctime(), d)
	def bcb(d):
		for (k, v) in d.items():
			cb("{} : {}".format(k, v))
	print assemble(path='/home/justin.payne', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN003966/CFSAN003966_01/CFSAN003966_S6_L001_R1_001.fastq', 
			 reads2='/shared/gn2/CFSANgenomes/CFSAN003966/CFSAN003966_01/CFSAN003966_S6_L001_R2_001.fastq', 
			 accession='CFSAN003966_01',  
			 insert_size=500,
			 k_value=32,
			 callback=cb,
			 update_callback=bcb,
			 debug='-debug' in sys.argv)