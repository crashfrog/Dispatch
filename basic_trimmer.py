import subprocess
import assembly_resources as tempfile
import shutil
import os

"""

	Basic first-16-bases trimming for MiSeq FASTQ's.

"""

def trim(assembler=lambda d: None, reads1=None, reads2=None, trimmer_args="", callback=lambda s, t: None, **kwargs):
	"Trims the weird GC-ratio region from each read; this is the first 16 bases."
	if trimmer_args:
		trimmer_args = " " + trimmer_args #pad with a space
	
	temp_path = tempfile.mkdtemp()
	temp_reads1 = os.path.join(temp_path, os.path.basename(reads1).replace(".fastq", ".basic_trimming.fastq"))
	
	
	callback("running trimming ({0})".format(reads1), status='Running trimming')
	#print "trimming ", temp_reads1, temp_reads2
	subprocess.check_call("fastx_trimmer -Q 33 -t 16{0} -i {1} -o {2}".format(trimmer_args, reads1, temp_reads1), shell=True)
	kwargs['reads1'] = temp_reads1
	
	if reads2:
		#print reads2
		callback("trimming ({0})".format(reads2), status='Running trimming')
		temp_reads2 = os.path.join(temp_path, os.path.basename(reads2).replace(".fastq", ".basic_trimming.fastq"))
		subprocess.check_call("fastx_trimmer -Q 33 -t 16{0} -i {1} -o {2}".format(trimmer_args, reads2, temp_reads2), shell=True)
		kwargs['reads2'] = temp_reads2
	
	
	callback("trimming done ({0})".format(temp_path))
	
	
	
	kwargs['callback'] = callback
	
	try:
		d = assembler.assemble(**kwargs)
	
	finally:
		try:
			callback("Cleaning up trimmed reads: {0}".format(temp_path))
			shutil.rmtree(temp_path)
		except Exception:
			if kwargs['debug']:
				import traceback, sys
				traceback.print_exc(sys.stdout)
	
	return d