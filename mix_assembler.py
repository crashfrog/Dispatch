#!/usr/bin/env python

import subprocess
import tempfile
import shutil
import os
import re

import fasta_statter


"""
	
	Mix, the assembly-aggregating microbial genome finishing tool.
	
	See: "Finishing bacterial genome assemblies with Mix"
	Soueidan et al. BMC Bioinformatics 2013, 14(Suppl 15):S16
	http://www.biomedcentral.com/1471-2105/14/S15/S16

"""
description = "Mix, the assembly-aggregating microbial genome finishing tool. Runs everything and then merges the assemblies via contig extension."

core_load = 8 #number of cores this assembler will max out

dont_use = ("Supersembler", "Velvet", "WORST", "Mix")

supports = ('MiSeq','IonTorrent','PacBio','454')

path = "/home/justin.payne/MIX/bin/"

def assemble(data_type, ref_file=None, assembler_dict={}, callback=lambda s: None, update_callback=lambda d: None, debug=False, **kwargs):
	d = {'assembler':"Mix",
		 'assembly_version':'Mix 1.0',
		 'average_coverage':'',
		 'num_contigs':'',
		 'n50':'',
		 'num_bases':'',
		 'fasta_file':'',
		 'lib_insert_length':'Not determined',
		 'matched':''
		 }
		 
	curr_dir = os.getcwd()
	try:
		working_dir = tempfile.mkdtemp()
		os.chdir(working_dir)
		results = list()
		for assembler_name, assembler in assembler_dict.items():
			if assembler_name not in dont_use:
				def callback_wrapper(s): #closure to wrap callback method, if any, from assembly_dispatch
					callback("{}: {}".format(assembler_name, s))
				if data_type in assembler.supports:
					callback("Looking for {}...".format(os.path.join(kwargs['path'], "{}.{}.fasta".format(kwargs['accession'], assembler_name))))
					if not os.path.exists(os.path.join(kwargs['path'], "{}.{}.fasta".format(kwargs['accession'], assembler_name))):
						try:
							args = {}
							args.update(kwargs)
							args['path'] = working_dir
							args['debug'] = debug
							r = assembler.assemble(data_type=data_type, fasta_file_name="{}.{}.fasta".format(kwargs['accession'], assembler_name), callback=callback_wrapper, debug=debug, **args)
							r['assembler'] = assembler_name
							results.append(r)
						except Exception:
							import traceback
							import sys
							traceback.print_exc(sys.stdout)
					else:
						shutil.copyfile(os.path.join(kwargs['path'], "{}.{}.fasta".format(kwargs['accession'], assembler_name)), os.path.join(working_dir, "{}.{}.fasta".format(kwargs['accession'], assembler_name)))
						results.append({'assembler':assembler_name, 'fasta_file':"{}.{}.fasta".format(kwargs['accession'], assembler_name)})
					
		if (not len(results)) or (len(results) < 2 and not debug):
			raise ValueError("No assembly returned results.")
		callback("Mix: running preprocessing") 
		print subprocess.check_output("{}preprocessing.py -o {}/{}.nucmer.fasta {}".format(
			path,
			working_dir,
			kwargs['accession'],
			" ".join([r['fasta_file'] for r in results])
		)
		, shell=True)
		
		
		callback("Mix: running NUCmer")
		subprocess.check_call("""nucmer --maxmatch --banded -c 30 -l 30 -p "alignments" {working_dir}/{accession}.nucmer.fasta {working_dir}/{accession}.nucmer.fasta""".format(working_dir=working_dir, **kwargs), shell=True)
		callback("Mix: running coords")
		subprocess.check_call("show-coords -rcl alignments.delta > alignments.coords", shell=True)
		callback("Mix: running mix")
		if not os.path.exists("Mix"):
			os.mkdir("Mix")
		subprocess.check_call("{}Mix.py -a alignments.coords -c {}.nucmer.fasta -o Mix/ -A 500 -C 0".format(path, kwargs['accession']), shell=True)
		shutil.copyfile("Mix/Mix_results_A500_C0/Mix_assembly.fasta", "{path}/{accession}.mix.fasta".format(**kwargs))
		
		d.update(fasta_statter.stat_velvet("{path}/{accession}.mix.fasta".format(**kwargs)))
		
		results.append(d)
		
		quast_results = fasta_statter.quast_compare(kwargs['path'], results, callback=callback, update_callback=update_callback, gi=ref_file, debug=debug)
		
	
	finally:
		os.chdir(curr_dir)
		if not debug:
			callback("Cleaning temp dir {}...".format(working_dir))
			shutil.rmtree(working_dir)
		
	return d
	
	
if __name__ == "__main__":
	#debug
	import datetime
	def cb(d):
		print "[{}] {}".format(datetime.datetime.today().ctime(), d)
	def bcb(d):
		for (k, v) in d.items():
			cb("{} : {}".format(k, v))
	d = dict()
	import assembly_dispatch
	d.update(assembly_dispatch.assembler_dict)
	del d['CLC']
	del d['Celera']
	print assemble(path='/home/justin.payne', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN003966/CFSAN003966_01/CFSAN003966_S6_L001_R1_001.fastq', 
			 reads2='/shared/gn2/CFSANgenomes/CFSAN003966/CFSAN003966_01/CFSAN003966_S6_L001_R2_001.fastq', 
			 accession='CFSAN003966_01',  
			 callback=cb,
			 update_callback=bcb,
			 k_value=144,
			 insert_size=500,
			 debug=True,
			 data_type = "MiSeq",
			 assembler_dict=d)
# 	print assemble(path='/home/justin.payne', 
# 			 reads1='/shared/gn2/CFSANgenomes/CFSAN001812/CFSAN001812_01/R_2012_04_12_16_25_14_user_IT1-33-BW3_Auto_IT1-33-BW3_38.fastq', 
# 			 accession='CFSAN001812_01', 
# 			 callback=cb,
# 			 update_callback=bcb,
# 			 k_value=144,
# 			 insert_size=500,
# 			 data_type='IonTorrent',
# 			 debug=True,
# 			 assembler_dict=d)