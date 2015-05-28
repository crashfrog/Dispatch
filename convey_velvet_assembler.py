import subprocess
import re
import datetime
import shutil
import fasta_statter
import sys


if "check_output" not in dir( subprocess ): # duck punch it in!
    def f(*popenargs, **kwargs):
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd)
        return output
    subprocess.check_output = f

ins_expression = r"Paired-end library 3 has length: (?P<ins_length>\d*)"
exp_cov_expression = r"\[\d*\.\d*\] Estimated Coverage = (?P<coverage>\d*\.\d*)"
exp_cut_expression = r"\[\d*\.\d*\] Estimated Coverage cutoff = (?P<coverage_cutoff>\d*\.\d*)"
parse_expression = r"""Final graph has (?P<num_nodes>\d*) nodes and n50 of (?P<n50>\d*), max (?P<max>\d*), total (?P<num_bases>\d*), using (?P<num_reads>\d*)/(?P<total_reads>\d*) reads"""

description = "Convey's accellerated version of the canonical de Brujin graph assembler. Parameterized for paired-end or single-end sequencing."

core_load = 1 #number of CPU cores this assembler will max out

supports = ('MiSeq','IonTorrent')

def assemble(accession, path, reads1, reads2=None, insert_size=500, minimum_contig_length=200, k_value=51, exp_coverage='auto', debug=True, callback=lambda s: None, fasta_file_name=None, **kwargs):

	#raise NotImplementedError()

	start_time = datetime.datetime.now()
	
	#get velvet version info
	#version = re.search(r"Version (\d*\.\d*\.\d*)", subprocess.check_output("velveth --version", shell=True)).groups()[0]
	
	#callback(subprocess.check_output("cnyinfo", shell=True))
	[callback(l) for l in subprocess.check_output("cnyinfo", shell=True).split("\n")]

	if reads2:
		velveth = "cnygc {0}/{1} {2} -fastq -shortPaired -separate {3} {4}"
		velveth_expression = velveth.format(path, accession, k_value, reads1, reads2)
	else:
		velveth = "cnygc {0}/{1} {2} -fastq -short {3}"
		velveth_expression = velveth.format(path, accession, k_value, reads1)
	
	velvetg = "velvetg {0}/{1} -ins_length {2} -min_contig_lgth {3} -exp_cov {4} -cov_cutoff auto -very_clean yes -scaffolding no"
	velvetg_expression = velvetg.format(path, accession, insert_size, minimum_contig_length, exp_coverage)
	
	#print "assembling ", reads1, reads2
	callback('running cnygc')
	if not debug:
		import os
		buffer = open(os.devnull, 'w')
	else:
		buffer = sys.stdout
		
	buffer.write(subprocess.check_output(velveth_expression, shell=True))
	
	callback('running velvetg')

	result = subprocess.check_output(velvetg_expression, shell=True)
	buffer.write(result)

	
	work_time = datetime.datetime.now() - start_time
	#capture results and look for assembly characteristics
	if not fasta_file_name:
		fasta_file_name = "{0}.fasta".format(accession)
	try:
		m = re.search(parse_expression, result)
		if m:
			d = m.groupdict()
			d['assembly_version'] = 'Velvet v. {0}'.format(version)
			
			if "EMPTY GRAPH" not in result:

				shutil.move("{}/{}/contigs.fa".format(path, accession), "{0}/{1}".format(path, fasta_file_name))
		
				c = re.search(ins_expression, result)
				if c:
					cd = c.groupdict()
					d['lib_insert_length'] = cd['ins_length']
				else:
					d['lib_insert_length'] = 'Not a number'
			

				d.update(fasta_statter.stat_velvet("{0}/{1}".format(path, fasta_file_name), k_value))
	
				d['fasta_file'] = fasta_file_name
				
			else:
				raise ValueError("Assembly produced empty graph.")
	finally:
		with open("{0}/{1}.stats".format(path, accession), 'w') as stats:
			stats.write("""Automated Velvet {ver} assembly
completed {date} in {time}
Called:
{velveth}
{velvetg}
{output}
""".format(ver=version,
			   date=datetime.datetime.today().isoformat(),
			   time=str(work_time),
			   velveth=velveth_expression,
			   velvetg=velvetg_expression,
			   output=result))
		   
	
		   
	
	
	return d
	
if __name__ == "__main__":
	#debug
	import datetime
	def cb(d):
		print "[{0}] {1}".format(datetime.datetime.today().ctime(), d)
	def bcb(d):
		for (k, v) in d.items():
			cb("{0} : {1}".format(k, v))
	print assemble(path='/home/justin.payne', 
			 reads1='/home/justin.payne/mounts/dna/gnome2/CFSANgenomes/CFSAN003905/CFSAN003905_01/CFSAN003905_S12_L001_R1_001.fastq', 
			 reads2='/home/justin.payne/mounts/dna/gnome2/CFSANgenomes/CFSAN003905/CFSAN003905_01/CFSAN003905_S12_L001_R2_001.fastq', 
			 accession='CFSAN003905_01', 
			 callback=cb,
			 update_callback=bcb,
			 k_value=31)
	print assemble(path='/home/justin.payne', 
			 reads1='/home/justin.payne/mounts/dna/gnome2/CFSANgenomes/CFSAN001812/CFSAN001812_01/R_2012_04_12_16_25_14_user_IT1-33-BW3_Auto_IT1-33-BW3_38.fastq', 
			 accession='CFSAN001812_01', 
			 callback=cb,
			 update_callback=bcb,
			 )

