import subprocess
import tempfile
import shutil
import os
from Bio import SeqIO

"""

	"Diginorm" trimming: eliminate read pairs if neither pair contributes unique k-mers to a kmer-set. Uses same k value
	as underlying assembler.

"""

def substrings(s, k):
	"Return a list of all of the k-mer substrings in a string."
	return [s[i:i+k] for i in range(0, len(s) - k)]
		

def trim(assembler=lambda d: None, reads1=None, reads2=None, trimmer_args="", k_value="61", callback=lambda s, t: None, **kwargs):
	"Trims the weird GC-ratio region from each read; this is the first 16 bases."
	if trimmer_args:
		trimmer_args = " " + trimmer_args #pad with a space
	
	try:
		temp_path = tempfile.mkdtemp()
		temp_reads1 = os.path.join(temp_path, os.path.basename(reads1).replace(".fastq", ".diginorm_trimmed.fastq"))
	
	
		callback("trimming ({})".format(reads1), status='Running trimming')
	
		kwargs['reads1'] = temp_reads1
		kwargs['k_value'] = k_value = int(k_value)
	
		if reads2:
			#print reads2
			callback("trimming ({})".format(reads2), status='Running trimming')
			temp_reads2 = os.path.join(temp_path, os.path.basename(reads2).replace(".fastq", ".diginorm_trimmed.fastq"))
			kwargs['reads2'] = temp_reads2
		
			kmer_set = {} #it may be faster to add to a dict than to a set
		
			with open(reads1, 'r') as rin_1:
				with open(reads2, 'r') as rin_2:
					with open(temp_reads1, 'w') as rout_1:
						with open(temp_reads2, 'w') as rout_2:
							#diginorm trimming - reject a read if it adds no distinct k-mers to the k-mer set
							for contig1, contig2 in zip(SeqIO.parse(rin_1, 'fastq'), SeqIO.parse(rin_2, 'fastq')):
								start_size = len(kmer_set)
								for kmer in substrings(str(contig1.seq), k_value) + substrings(str(contig2.seq), k_value):
									kmer_set[kmer] = None
								if len(kmer_set) > start_size:
									rout_1.write(contig1.format('fastq'))
									rout_2.write(contig2.format('fastq'))
		else:
			kmer_set = {} #it may be faster to add to a dict than to a set
			with open(reads1, 'r') as rin_1:
				with open(temp_reads1, 'w') as rout_1:
					for contig1 in SeqIO.parse(rin_1, 'fastq'):
						start_size = len(kmer_set)
						for kmer in substrings(str(contig1.seq), k_value):
							kmer_set[kmer] = None
						if len(kmer_set) > start_size:
							rout_1.write(contig1.format('fastq'))
	
	
		callback("trimming done ({})".format(temp_path))
	
	
	
		kwargs['callback'] = callback
	
		d = assembler.assemble(**kwargs)
	
	finally:
		try:
			callback("Cleaning up trimmed reads: {}".format(temp_path))
			shutil.rmtree(temp_path)
		except Exception:
			pass
	
	return d