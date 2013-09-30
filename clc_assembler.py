#!/usr/bin/env python

import subprocess
import re
import datetime
import shutil
import os

"""
	CLC Server CLI binding for Pipeline: Assembly Dispatch.
	Not meant to run on its own.

"""

description = "CLC Bio's proprietary de Brujin assembler; runs on CLC Genomics Server."

core_load = 8 #number of CPU cores this assembler will use (doesn't really use any, but lets not go crazy)

supports = ('MiSeq', 'IonTorrent')

output_expression = r'File: (?P<url>[\w\W]*?)\n//\n'

file_id_expression = r"ClcUrl: (?P<url>[\w\W]*?)\n//\n"

base_parameters = {'user':'justin.payne',
				   'token':'BAAAAAAAAAAAAAP583c10bfdbd326ba-32830091-14113965181--8000'  #generated using clcserverkeystore
				   } 

def assemble(path, data_type='MiSeq', callback=lambda s: None, **kwargs):

	os.chdir('/shared') #a lot of shell manipulation in this, so it's important to start at a known location

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
	
	subprocess.check_output("clcserver -S xserve13.fda.gov -U {user} -W {token} -A mkdir -t clc://server/CLCData/AutomatedAssembly/Reads/ -n {accession}".format(**base_parameters), shell=True)
	subprocess.check_output("clcserver -S xserve13.fda.gov -U {user} -W {token} -A mkdir -t clc://server/CLCData/AutomatedAssembly/Assemblies/ -n {accession}".format(**base_parameters), shell=True)
	
	
	
	callback('running import')
	
	
	#print os.getcwd()
	
	try:
		if 'MiSeq' in data_type:
			results = subprocess.check_output("clcserver -S xserve13.fda.gov -U {user} -W {token} -A ngs_import_illumina --paired-reads true --quality-score ILLUMINA_PHRED "
									"-d clc://server/CLCData/AutomatedAssembly/Reads/{accession}/ "
									"-f {reads1} -f {reads2}".format(**base_parameters), shell=True)
		elif 'IonTorrent' in data_type:
			results = subprocess.check_output("clcserver -S xserve13.fda.gov -U {user} -W {token} -A ngs_import_iontorrent "
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
		
								
		results = subprocess.check_output("clcserver -S xserve13.fda.gov -U {user} -W {token} -A denovo_assembly "
								"--perform-scaffolding false "
								"-i {url} "
								"-d clc://server/CLCData/AutomatedAssembly/Assemblies/{accession}/ ".format(**base_parameters), shell=True)
								
		m = re.findall(file_id_expression, results)
		if m:
			base_parameters['report'] = m[0]
			base_parameters['assembly'] = m[1]
			
			os.chdir(base_parameters['asm'])
			results = subprocess.check_output("clcserver -S xserve13.fda.gov -U {user} -W {token} -A export "
									"-s {report} -d {asm} -e excel_1997_2007".format(**base_parameters), shell=True)
									
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
		
									
		results = subprocess.check_output("clcserver -S xserve13.fda.gov -U {user} -W {token} -A export "
									"-s {assembly} -e fasta -d {asm}".format(**base_parameters), shell=True)
									
		#print results
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
			subprocess.call("clcserver -S xserve13.fda.gov -U {user} -W {token} -A rm -r -t clc://server/CLCData/AutomatedAssembly/Reads/{accession}/ > /dev/null".format(**base_parameters), shell=True)
			subprocess.call("clcserver -S xserve13.fda.gov -U {user} -W {token} -A rm -r -t clc://server/CLCData/AutomatedAssembly/Assemblies/{accession}/ > /dev/null".format(**base_parameters), shell=True)

		except subprocess.CalledProcessError as e:
			print e
			print base_parameters
			print d
	os.chdir('/shared')
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
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001812/CFSAN001812_01/R_2012_04_12_16_25_14_user_IT1-33-BW3_Auto_IT1-33-BW3_38.fastq', 
			 accession='CFSAN001812_01', 
			 callback=cb,
			 update_callback=bcb,
			 data_type='IonTorrent')
#clc_assembler.assemble(path='gn2/CFSANgenomes/CFSAN004196/asm/', reads1='gn2/CFSANgenomes/CFSAN004196/CFSAN004196_01/CFSAN004196_S9_L001_R1_001.fastq', reads2='gn2/CFSANgenomes/CFSAN004196/CFSAN004196_01/CFSAN004196_S9_L001_R2_001.fastq', accession='CFSAN004196_01', callback=lambda s: print s)