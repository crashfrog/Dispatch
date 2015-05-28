import subprocess
import tempfile
import shutil
import os.path as p
import re
import sys

import xmlrpclib
from Bio import SeqIO

import fasta_statter

"""
	
	Uses three strategies to pick reference:
		1) Use reference specified by job.
		2) Look it up based on species.
		3) Do a quick Velvet and find best match to first 1500-bp contig.

"""
description = "Auto-reference-picking assembler."

core_load = 8 #number of cores this assembler will max out

supports = ('MiSeq','IonTorrent','PacBio','454')

def assemble(accession, ref_file=None, callback=lambda s: None, update_callback=lambda d: None, debug=False, **kwargs)
	d = {'assembly_version':'',
		 'average_coverage':'',
		 'num_contigs':'',
		 'n50':'',
		 'num_bases':'',
		 'fasta_file':'',
		 'lib_insert_length':'Not determined',
		 'matched':''
		 }
	working_dir = tempfile.mkdtemp()
	try:
		#get version
		ver = re.search("Version: (\S+)\n", subprocess.check_output("bwa7", shell=True))
		if ver:
			d['assembly_version'] = "BWA v. {}".format(ver.group(1))
			
		if ref_file:
			urllib.urlretrieve('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text'.format(ref_file),
							   os.path.join(working_dir, "reference.fasta")
		else:
			try:
				server = xmlrpclib.ServerProxy("http://cfe1019692:8080")
				data = server.retrieve_cfsan(accession)
				name = "{Genus} {Species} subsp. {Subspecies} serovar {Serovar} str. {StrainName}".format(**data).replace("serovar str.", "str.").replace("subsp. str.", "str.")
				#query Entrez for first-match reference
				
				update_callback({"ref_file":something, "ref_url":something})
			except Exception:					   
				working_dir = tempfile.mkdtemp()
				import velvet_assembler_notrimmer
				r = velvet_assembler_notrimmer.assemble(accession=accession,
														path=working_dir,
														callback=callback,
														update_callback=update_callback,
														k_value=63,
														debug=debug,
														fasta_file_name="{}.velvet.fasta".format(accession)
														**kwargs)
				ref_file = fasta_statter.find_closest_ref(os.path.join(working_dir, r['fasta_file']), callback=lambda s: callback("BWA: {}".format(s)), update_callback=update_callback)
				urllib.urlretrieve('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text'.format(ref_file),
								   os.path.join(working_dir, "reference.fasta")
		
		out = open("/dev/null", 'w')
		if debug:
			out = sys.stdout
		callback("Running: indexing reference")
		out.write(subprocess.check_output("bwa7 index -a is {}/reference.fasta".format(working_dir), shell=True))
		callback("Running: mapping reads to reference")
		contigs = SeqIO.parse(subprocess.check_output("bwa7 mem {}/reference.fasta {reads1} {reads2}
		#assemble
	
	d.update(fasta_statter.stat_fasta(...))
	
	finally:
		try:
			if not debug:
				shutil.rmtree(working_dir)
		except:
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