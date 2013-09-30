#!/usr/bin/env python
import shutil
import re
import os
import datetime
import subprocess
import tempfile
import sys

import fasta_statter

"""

	Assembly Dispatch binding for the hoary old Celera Assembler. A good template for
	future assembly bindings.

"""

description = "WGS/Celera assembler, the overlap-layout-consensus workhorse for long reads."

core_load = 1 #number of CPU cores this assembler will max out

supports = ('MiSeq','IonTorrent','PacBio','454')

def assemble(accession, path, reads1,  insert_size=400, reads2=None, data_type='MiSeq', callback=lambda s: None, fasta_file_name=None, debug=True, **kwargs):
	working_path = tempfile.mkdtemp()
	#raise ValueError("Celera assembler not currently supported.")
	try:
		d = {'assembly_version':'Celera Assembler v. 7',
			 'average_coverage':'',
			 'num_contigs':'',
			 'n50':'',
			 'num_bases':'',
			 'fasta_file':'',
			 'lib_insert_length':'',
			 'matched':''
			 }
		if not fasta_file_name:
			fasta_file_name = "{}.celera.fastq".format(accession)
		d['fasta_file'] = fasta_file_name
		
		if debug:
			buffer = sys.stdout
		else:
			buffer = open(os.devnull, 'w')
	
		asm_path = os.path.join(path, '../asm')
		

		frg = os.path.join(working_path, "{}.frg".format(accession))
		subprocess.call('export PATH="$PATH:/usr/local/bin/Linux-amd64/bin"', shell=True)
	
		#interleave reads and produce FRG file
		callback('running fastqToCA')
		if 'MiSeq' in data_type:
			buffer.write(subprocess.check_output("fastqToCA -insertsize {} 20 -libraryname {} -technology illumina -type illumina -mates {},{} > {}".format(insert_size, accession, reads1, reads2, frg), shell=True))
		elif 'PacBio' in data_type:
			raise ValueError("PacBio data not yet supported.")
		else: #454 or Iontorrent
			buffer.write(subprocess.check_output("fastqToCA -insertsize {} 20 -libraryname {} -technology 454 -reads {} > {}".format(insert_size, accession, reads1, frg), shell=True))

		callback('running celera')
		results = subprocess.check_output("runCA -d {} -p {} gkpFixInsertSizes=0 {}".format(working_path, accession, frg), shell=True)
		buffer.write(results)
		
		m = re.search(r"TotalContigsInScaffolds\s*(\d*)", results)
		if m:
			d['num_contigs'] = m.group(1)
			
		m = re.search(r"TotalBasesInScaffolds\s*(\d*)", results)
		if m:
			d['num_bases'] = m.group(1)
			
		m = re.search(r"N50ContigBases\s*(\d*)", results)
		if m:
			d['n50'] = m.group(1)
			
		m = re.search(r"ContigsOnly\s*(\d*.\d*)", results)
		if m:
			d['average_coverage'] = str(int(round(float(m.group(1))))) + "X"
			
		buffer.write(subprocess.check_output("asmOutputFasta -C -D -S -p {}/{}.celera < {}/{}.asm".format(path, accession, working_path, accession), shell=True))
		for token in ("utg.qual", "utg.qv"):
			os.remove(os.path.join(path, "{}.celera.{}".format(accession, token)))
		shutil.copyfile("{}/{}.celera.utg.fasta".format(path, accession), "{}/{}.celera.fasta".format(path, accession))
		d.update(fasta_statter.stat_fasta("{}/{}.celera.fasta".format(path, accession)))
		
	finally:
		callback("Cleaning up {}...".format(working_path))
		shutil.rmtree(working_path)

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
			 reads1='/shared/gn2/CFSANgenomes/CFSAN003966/CFSAN003966_01/CFSAN003966_S6_L001_R1_001.fastq', 
			 #reads2='/shared/gn2/CFSANgenomes/CFSAN003966/CFSAN003966_01/CFSAN003966_S6_L001_R2_001.fastq', 
			 accession='CFSAN003966_01',
			 data_type='IonTorrent',
			 insert_size=500,
			 callback=cb,
			 update_callback=bcb,
			 debug=True)