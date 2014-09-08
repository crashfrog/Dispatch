#!/usr/bin/env python

import subprocess
import re
import datetime
import shutil
import os

"""
	CLC Server CLI binding for Pipeline: Assembly Dispatch.
	Not meant to run on its own, but can.

"""

usage = """
Usage: clc_assembler [options] reads1.fastq [reads2.fastq] [-o /path/to/output]

Options
	-o	Output directory

"""

description = "CLC Assembly Cell proprietary de Brujin assembler; runs on CLC Genomics Server."

core_load = 8 #number of CPU cores this assembler will use (doesn't really use any, but lets not go crazy)

supports = ('MiSeq', 'IonTorrent', 'MiSeq_MP')

output_expression = r'File: (?P<url>[\w\W]*?)\n//\n'

file_id_expression = r"ClcUrl: (?P<url>[\w\W]*?)\n//\n"

base_parameters = {'user':'justin.payne',
				   'token':'BAAAAAAAAAAAAAP583c10bfdbd326ba--76de074e-1449e8a81d4--8000',  #generated using clcserverkeystore
				   'server':'xserve13.fda.gov'
				   } 

def assemble(path, data_type='MiSeq', callback=lambda s: None, **kwargs):


	d = {'assembly_version':'CLC Genomics Server v. 5.0.1',
		 'average_coverage':'',
		 'num_contigs':'',
		 'n50':'',
		 'num_bases':'',
		 'fasta_file':'',
		 'lib_insert_length':'not determined',
		 'matched':''
		 }
		 
	
	
	base_parameters.update(kwargs)
	base_parameters['asm'] = path
	
	callback('running directory creation')
	
	subprocess.check_output("clcserver -S {server} -U {user} -W {token} -A mkdir -t clc://server/CLCData/AutomatedAssembly/Reads/ -n {accession}".format(**base_parameters), shell=True)
	subprocess.check_output("clcserver -S {server} -U {user} -W {token} -A mkdir -t clc://server/CLCData/AutomatedAssembly/Assemblies/ -n {accession}".format(**base_parameters), shell=True)
	
	
	
	callback('running import')
	
	
	#print os.getcwd()
	
	try:
		if 'MiSeq' in data_type:
			results = subprocess.check_output("clcserver -S {server} -U {user} -W {token} -A ngs_import_illumina --paired-reads true --quality-score ILLUMINA_PHRED "
									"-d clc://server/CLCData/AutomatedAssembly/Reads/{accession}/ "
									"-f {reads1} -f {reads2}".format(**base_parameters), shell=True)
		elif 'IonTorrent' in data_type:
			results = subprocess.check_output("clcserver -S {server} -U {user} -W {token} -A ngs_import_iontorrent "
									"-d clc://server/CLCData/AutomatedAssembly/Reads/{accession}/ "
									"-f {reads1}".format(**base_parameters), shell=True)
		else:
			raise ValueError('Data type "{}" not currently supported.'.format(data_type))
	
							
		#print results
							
		m = re.search(file_id_expression, results)
		if m:
			base_parameters['url'] = m.group('url')
		else:
			print "Output:", results
			raise ValueError()
	
			
		callback('running CLC assembler')
		
								
		results = subprocess.check_output("clcserver -S {server} -U {user} -W {token} -A denovo_assembly "
								"--perform-scaffolding false "
								"--create-report-1 true --create-report-2 true --create-report-3 true "
								"--map-reads-to-contigs MAP_BACK -i {url} "
								"-d clc://server/CLCData/AutomatedAssembly/Assemblies/{accession}/ ".format(**base_parameters), shell=True)
								
		m = re.findall(file_id_expression, results)
		if m:
			base_parameters['report'] = m[0]
			base_parameters['assembly'] = m[1]
			
			os.chdir(base_parameters['asm'])
			
# 			results = subprocess.check_output("clcserver -S {server} -U {user} -W {token} "
# 											  "-A detailed_mapping_report"
			
			callback("parsing report")
			results = subprocess.check_output("clcserver -S {server} -U {user} -W {token} -A export "
									"-i {report} -d {asm} -e excel_1997_2007".format(**base_parameters), shell=True)
									
			m = re.search(output_expression, results)
			if m:
				d['report'] = m.group('url')
				import xlrd
				x = xlrd.open_workbook(os.path.join(base_parameters['asm'], m.group('url')))
				for key, value, extra in zip(x.sheets()[0].col_values(0), x.sheets()[0].col_values(1), x.sheets()[0].col_values(3)):
					if 'N50' in str(key):
						d['n50'] = str(int(value))
					elif 'Count' in str(key):
						d['num_contigs'] = str(int(value))
					elif 'Total' in str(key):
						d['num_bases'] = int(value)
					elif 'Matched' in str(key):
						d['matched'] = int(extra)
				try:
					d['average_coverage'] = '{}X'.format(str(int(round(d['matched'] / d['num_bases']))))
				except (ZeroDivisionError, ValueError, TypeError):
					d['average_coverage'] = 'not a number'

						
		else:
			print "Output:", results
			raise ValueError('Assembly returned no results.')
				
		callback('running assembly retrieval')
		
									
		results = subprocess.check_output("clcserver -S {server} -U {user} -W {token} -A export "
									"-i {assembly} -e fasta -d {asm}".format(**base_parameters), shell=True)
									
		m = re.search(output_expression, results)
		if m:
			fasta = m.group('url')
			if 'fasta_file_name' in base_parameters:
				d['fasta_file'] = base_parameters['fasta_file_name']
			else:
				d['fasta_file'] = "{accession}.clc.fasta".format(**base_parameters)
			shutil.move(fasta, os.path.join(base_parameters['asm'], d['fasta_file'])) #rename .fa file

		#print d

			
		
	finally:
		try:
			subprocess.call("clcserver -S {server} -U {user} -W {token} -A rm -r -t clc://server/CLCData/AutomatedAssembly/Reads/{accession}/ > /dev/null".format(**base_parameters), shell=True)
			subprocess.call("clcserver -S {server} -U {user} -W {token} -A rm -r -t clc://server/CLCData/AutomatedAssembly/Assemblies/{accession}/ > /dev/null".format(**base_parameters), shell=True)

		except subprocess.CalledProcessError as e:
			print e
			print base_parameters
			print d

	return d
			
if __name__ == "__main__":
	
	import datetime
	import sys
	def cb(d):
		print "[{}] {}".format(datetime.datetime.today().ctime(), d)
	def bcb(d):
		for (k, v) in d.items():
			cb("{} : {}".format(k, v))
	if len(sys.argv) == 1:
		print usage
		quit()
	else:
		if '-debug' in sys.argv:
			print assemble(path='/home/justin.payne', 
					 reads1='/shared/gn2/CFSANgenomes/CFSAN001812/CFSAN001812_01/R_2012_04_12_16_25_14_user_IT1-33-BW3_Auto_IT1-33-BW3_38.fastq', 
					 accession='CFSAN001812_01', 
					 callback=cb,
					 update_callback=bcb,
					 data_type='IonTorrent')
			quit()
		path = os.getcwd()
		if '-o' in sys.argv:
			#output directory
			path = sys.argv.pop(sys.argv.index('-o') + 1)
			sys.argv.remove(sys.argv.index('-o'))
		reads2 = None
		data_type = 'IonTorrent'
		if len(sys.argv) > 2:
			reads2 = os.path.abspath(sys.argv[2])
			data_type='MiSeq'
		reads1 = os.path.abspath(sys.argv[1])
		bcb(assemble(path=path,
					 accession=reads1.split("/")[-1].split(".")[0].split("_")[0],
					 reads1=reads1,
					 reads2=reads2,
					 callback=cb,
					 data_type=data_type))