#!/usr/bin/env python

import subprocess
import tempfile
import shutil
import os

"""

	Basic first-16-bases trimming for MiSeq FASTQ's, plus kmer analysis by KmerGenie.
	
	Bioinformatics. 2014 Jan 1;30(1):31-7. doi: 10.1093/bioinformatics/btt310. Epub 2013 Jun 3.
	Informed and automated k-mer size selection for genome assembly.
	Chikhi R, Medvedev P.

"""

def trim(assembler=lambda d: None, reads1=None, reads2=None, trimmer_args="", callback=lambda s, t: None, update_callback=lambda d: None, **kwargs):
	"Trims the weird GC-ratio region from each read; this is the first 16 bases. Then use KmerGenie to choose the right kmer size"
	if trimmer_args:
		trimmer_args = " " + trimmer_args #pad with a space
	
	temp_path = tempfile.mkdtemp()
	callback("Created temp directory {}".format(temp_path))
	temp_reads1 = os.path.join(temp_path, os.path.basename(reads1).replace(".fastq", ".basic_trimming.fastq"))
	
	
	callback("running trimming ({})".format(reads1), status='Running trimming')
	#print "trimming ", temp_reads1, temp_reads2
	subprocess.check_call("fastx_trimmer -Q 33 -t 16{} -i {} -o {}".format(trimmer_args, reads1, temp_reads1), shell=True)
	kwargs['reads1'] = temp_reads1
	
	if reads2:
		print reads2
		callback("trimming ({})".format(reads2), status='Running trimming')
		temp_reads2 = os.path.join(temp_path, os.path.basename(reads2).replace(".fastq", ".basic_trimming.fastq"))
		subprocess.check_call("fastx_trimmer -Q 33 -t 16{} -i {} -o {}".format(trimmer_args, reads2, temp_reads2), shell=True)
		kwargs['reads2'] = temp_reads2
	
	
	callback("trimming done ({})".format(temp_path))
	callback("KmerGenie starting...")
	callback("""

	Bioinformatics. 2014 Jan 1;30(1):31-7. doi: 10.1093/bioinformatics/btt310. Epub 2013 Jun 3.
	Informed and automated k-mer size selection for genome assembly.
	Chikhi R, Medvedev P.
	
""")
	with tempfile.NamedTemporaryFile('w') as reads_list:
		for reads in ('reads1', 'reads2'):
			if kwargs.get(reads):
				reads_list.write(kwargs.get(reads) + '\n')
		reads_list.flush()
		
		if kwargs['debug']:
			print 'Reads parameter file:'
			print open(reads_list.name, 'r').read()

		output = subprocess.check_output("/usr/local/bin/kmergenie {} -k 251 -l 13 -o {}/kmergenie".format(reads_list.name, temp_path), shell=True)
		kmer = output.split('\n')[-2].split(' ')[-1]
		if kwargs['debug']:
			print output
	
		update_callback({'k_value':kmer, 'trimmer_args':'K-value selection by KmerGenie', 'trimmer':'kmergenie_trimmer'})
		
		try:
			shutil.copyfile(os.path.join(temp_path, 'kmergenie_report.html'), os.path.join(kwargs['path'], 'kmergenie_report.html'))
		except:
			if kwargs['debug']:
				import traceback
				traceback.print_exc()
	
		kwargs['callback'] = callback
		kwargs['update_callback'] = update_callback
	
	try:
		d = assembler.assemble(**kwargs)
	
	finally:
		try:
			callback("Cleaning up trimmed reads: {}".format(temp_path))
			if kwargs['debug']:
				raw_input("Debug hold; return to continue")
			shutil.rmtree(temp_path)
		except Exception:
			if kwargs['debug']:
				import traceback, sys
				traceback.print_exc(sys.stdout)
	
	return d
	

if __name__ == '__main__':
	#debug
	import datetime
	import worst_assembler
	def cb(d, *a, **v):
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