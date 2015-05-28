#!/usr/bin/env python

from Bio.SeqIO import *
import tempfile
import shutil
from os.path import join as j
import xml.etree.ElementTree as xml

def mean(l):
	return float(sum(l)) / float(len(l))

def trim(trimmer_args, reads1, reads2=None, callback=lambda s, **v: None, assembler=lambda d: None, debug=False, **kwargs):
	"Trim bases from low-quality side until a specified average quality is met."
	if debug:
		print trimmer_args
	params = xml.fromstring(trimmer_args).find('QUALITY_TRIMMER_ARGS')
	try:
		quality_cutoff = int(params.attrib['cutoff'])
	except KeyError:
		quality_cuttoff = 20
	callback('Starting trimming with cutoff q={}...'.format(quality_cutoff))
	
	temp_dir = tempfile.mkdtemp()
	callback('Running trimming in {}...'.format(temp_dir))
	try:
		if reads2:
		
			with open(reads1, 'r') as ri_1, open(reads2, 'r') as ri_2, open(j(temp_dir, 'reads1.fastq'), 'w') as ro_1, open(j(temp_dir, 'reads2.fastq'), 'w') as ro_2:
				num_reads = 0
				for (c1, c2) in zip(parse(ri_1, 'fastq'), parse(ri_2, 'fastq')):
					if mean(c1.letter_annotations['phred_quality']) > quality_cutoff and mean(c2.letter_annotations['phred_quality']) > quality_cutoff:
						ro_1.write(c1.format('fastq'))
						ro_2.write(c2.format('fastq'))
					num_reads += 1
					if num_reads % 100000 == 0:
						callback("Scanned {} reads...".format(num_reads))
		
			kwargs['reads1'] = j(temp_dir, 'reads1.fastq')
			kwargs['reads2'] = j(temp_dir, 'reads2.fastq')
			kwargs['callback'] = callback
		else:
			with open(reads1, 'r') as ri_1, open(j(temp_dir, 'reads1.fastq'), 'w') as ro_1:
				num_reads = 0
				for c1 in parse(ri_1, 'fastq'):
					if mean(c1.letter_annotations['phred_quality']) > quality_cutoff:
						ro_1.write(c1.format('fastq'))
					num_reads += 1
					if num_reads % 100000 == 0:
						callback("Scanned {} reads...".format(num_reads))
		
			kwargs['reads1'] = j(temp_dir, 'reads1.fastq')
			kwargs['callback'] = callback
	
		d = assembler.assemble(**kwargs)	
		
	finally:
		callback('Cleaning up {}...'.format(temp_dir))
		shutil.rmtree(temp_dir)
		
	return d
	
	
if __name__ == "__main__":
	import sys
	import datetime
	import abyss_assembler
	def cb(s):
		print "[{}]".format(datetime.datetime.today().ctime()), s
	def ucb(d):
		for (key, value) in d.items():
			print key, ":", value
	args = xml.Element('TRIMMER_ARGS')
	xml.SubElement(args, 'QUALITY_TRIMMER_ARGS', cutoff='23')
	ucb(trim(path='/home/justin.payne/', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001656/CFSAN001656_01/CFSAN001656_S8_L001_R1_001.fastq', 
			 reads2='/shared/gn2/CFSANgenomes/CFSAN001656/CFSAN001656_01/CFSAN001656_S8_L001_R2_001.fastq', 
			 accession='CFSAN001656_01',  
			 callback=cb,
			 debug='-debug' in sys.argv,
			 k_value='32',
			 trimmer_args=xml.tostring(args),
			 update_callback=ucb,
			 assembler=abyss_assembler))