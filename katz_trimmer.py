import subprocess
import tempfile
import shutil
import os

"""

	Lee Katz's (CDC) tool for quality-related trimming. This is a pretty good template
	for writing additional trimmers.

"""

def trim(assembler=lambda **k: dict(), reads1=None, reads2=None, trimmer_args="", callback=lambda s: None, **kwargs):
	"Lee Katz's MiSeq cleaning/trimming tool run_assembly_trimLowQualEnds. Trims by removing bases until GC skew is relatively constant."
	
	tempdir = tempfile.mkdtemp()
	try:
		if trimmer_args:
			trimmer_args += " "
		t1 = os.path.join(tempdir, "reads1_trimmed.fastq")
		t2 = os.path.join(tempdir, "reads2_trimmed.fastq")
		callback("Trimming {}".format(reads1))
		subprocess.check_call("run_assembly_trimLowQualEnds.pl {}{} > {}".format(trimmer_args, reads1, t1), shell=True)
		callback("Trimming {}".format(reads2))
		subprocess.check_call("run_assembly_trimLowQualEnds.pl {}{} > {}".format(trimmer_args, reads2, t2), shell=True)
	

		
	except subprocess.CalledProcessError as e:
		kwargs['reads1'] = reads1
		kwargs['reads2'] = reads2
		kwargs['callback'] = callback
	
		d = assembler(**kwargs)
	
	else:
		kwargs['reads1'] = t1
		kwargs['reads2'] = t2
		kwargs['callback'] = callback
	
		d = assembler(**kwargs)
	
	finally:
		shutil.rmtree(tempdir)
	
	return d
	
	
if __name__ == "__main__":
	def cb(s):
		print s
	print trim(reads1='/shared/gn2/CFSANgenomes/CFSAN004214/CFSAN004214_01/CFSAN004214_S10_L001_R1_001.fastq', 
			   reads2='/shared/gn2/CFSANgenomes/CFSAN004214/CFSAN004214_01/CFSAN004214_S10_L001_R2_001.fastq', 
			   trimmer_args="",
			   callback=cb)
			   
			   