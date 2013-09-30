#!/usr/bin/env python

import MySQLdb as mysql


import velvet_assembler_notrimmer

import velvet_optimize
import abyss_assembler
import clc_assembler
import wgs_assembler
import supersembler
import spades_assembler
#import masurca_assembler
import mira_assembler
import worst_assembler
#import idba_assembler

import basic_trimmer
import katz_trimmer
import quality_trimmer
import fasta_statter
import ten_each_end_trimmer

import datetime
from time import sleep
import subprocess
import shutil
import os
import sys

"""
	
	Assembly Dispatch Daemon
	May 1, 2013
	Justin Payne
	
	Assembly job dispatch daemon, single-threaded version. Interacts with Pipleline Job DB.
	Assembly bindings are written as submodules and encapsulate the shell commands necessary
	to operate the assemblers and capture statistics on an automated basis.
	
	This module handles the database interface and dispatches assembly jobs according to
	the "assembler_dict", below.

"""

def null_trimmer(assembler, **kwargs):
	"No trimming."
	return assembler.assemble(**kwargs)

assembler_dict = {'Velvet':velvet_assembler_notrimmer,
				  'Supersembler':supersembler,
				  'Celera':wgs_assembler,
				  'Velvet_optimize':velvet_optimize,
				  'ABySS':abyss_assembler,
				  'SPAdes':spades_assembler,
				  'CLC':clc_assembler,
				  #'MaSuRCA':masurca_assembler,
				  'Mira':mira_assembler,
				  #'IDBA-UD':idba_assembler,
				  'WORST':worst_assembler
				  }
				  
trimmer_dict = {'no_trimmer':null_trimmer, #trimmers can be in this module or others
				'basic_trimmer':basic_trimmer.trim,
				'katz_trimmer':katz_trimmer.trim, #currently unsupported
				'quality_trimmer':quality_trimmer.trim,
				'ten_each_end_trimmer':ten_each_end_trimmer.trim
				}



default_root = '/shared'


def query(s, commit_function=lambda l: True):
	"Generic method to abstract querying database; allows optional function to determine whether changes should be committed"
	items = []
	try:
		job_db = mysql.connect(host='xserve15.fda.gov', user='job_user', passwd='job_user', db='Jobs')
		jobc = job_db.cursor(mysql.cursors.DictCursor)
		jobc.execute(s)
		items = jobc.fetchall()
		if commit_function(items):
			job_db.commit()
		else:
			job_db.rollback()
		job_db.close()
	except mysql.DatabaseError as e:
		sys.__stdout__.write("[{}] {}:{}\n".format(datetime.datetime.today().ctime(), type(e), e))
	return items
				  
				  
def update_status_callback(id, message, status='Running',):
	sys.__stdout__.write("[{}] {}\n".format(datetime.datetime.today().ctime(), message)) #in case some assembler redirects print/stdout
	if len(message) > 80:
		query("UPDATE assemblies SET status='{status}', exception_note='{message}' WHERE id='{id}';".format(status=status.encode('string_escape'), message=message.encode('string_escape'), id=id))
	else:
		query("UPDATE assemblies SET status='{message}', exception_note='' WHERE id='{id}';".format(message=message.encode('string_escape'), id=id))


def assemble(job, debug=True):
	tempdir = None
	try:
		path = os.path.abspath(os.path.join(default_root, job['path_1'], '../asm'))
		if not os.path.exists(path):
			os.mkdir(path)
		try:
			exp_coverage = int(job['average_coverage'].replace('X', '').replace('x', ''))
		except ValueError:
			exp_coverage = 20
		assembler = assembler_dict.get(job['job_type'], velvet_assembler_notrimmer)
		trimmer = trimmer_dict.get(job['trimmer'], null_trimmer) #default - no trimming
		statter = getattr(fasta_statter, job['statter'])
		def update_callback(d):
			for (k, v) in d.items():
				print "[{}]\tmodule override: {}='{}'".format(datetime.datetime.today().ctime(), k, v)
				query("UPDATE assemblies SET {k}='{v}' WHERE id='{i}';".format(i=job['id'], k=k, v=str(v).encode('string_escape')))
		query("UPDATE assemblies SET status='running', exception_note='' WHERE id = '{}'".format(job['id']))
		
		if "gz" in job['file_1']:
			#unzip file to tempdir
			import tempfile
			import gzip
			tempdir = tempfile.mkdtemp()
			update_status_callback(job['id'], 'Unzipping reads in {}...'.format(tempdir))
			with gzip.open(os.path.join(default_root, job['path_1'], job['file_1'])) as r_in, open(os.path.join(tempdir, job['file_1'].replace(".gz", "")), 'w') as r_out:
				for l in r_in:
					r_out.write(l)
			if job['file_2']:
				with gzip.open(os.path.join(default_root, job['path_2'], job['file_2'])) as r_in, open(os.path.join(tempdir, job['file_2'].replace(".gz", "")), 'w') as r_out:
					for l in r_in:
						r_out.write(l)
			job['file_1'] = job['file_1'].replace(".gz", "")
			job['file_2'] = job['file_2'].replace(".gz", "")
			job['path_1'] = job['path_2'] = tempdir
		
		if 'sff' in job['file_1']:
			from Bio import SeqIO
			import tempfile
			tempdir = tempfile.mkdtemp()
			with open(os.path.join(default_root, job['path_1'], job['file_1']), 'rb') as r_in, open(os.path.join(tempdir, job['file_1'].replace(".sff", ".fastq")), 'w') as r_out:
				for c in SeqIO.parse(r_in, 'sff'):
					r_out.write(c.format('fastq'))
			job['file_1'] = job['file_1'].replace('.sff', '.fastq')
			job['path_1'] = tempdir
			
		reads1=os.path.join(default_root, job['path_1'], job['file_1'])
		reads2=None
		if job['file_2']:
			reads2=os.path.join(default_root, job['path_2'], job['file_2'])
		
		results = trimmer(assembler=assembler,
						  accession=job['accession'], path=path, 
						  reads1=reads1, 
						  reads2=reads2, 
						  insert_size=job['insert_size'] or 500, 
						  minimum_contig_length=job['min_contig_length'] or 200, 
						  k_value=job['k_value'] or 177, 
						  exp_coverage=exp_coverage or 'auto',
						  callback=lambda m, **s: update_status_callback(job['id'], m, **s),
						  assembler_dict=assembler_dict,
						  update_callback=update_callback,
						  statter=statter,
						  ref_file=job['ref_file'],
						  ref_url=job['ref_url'],
						  data_type=job['data_type'],
						  debug=debug,
						  trimmer_args=job['trimmer_args']
						  )
		try:
			job['job_type'] = 'Bowtie read mapping'
			results['lib_insert_length'] = fasta_statter.read_map(os.path.join(path, results['fasta_file']), 
																  os.path.join(default_root, job['path_1'], job['file_1']), 
																  os.path.join(default_root, job['path_2'], job['file_2']), 
																  callback=lambda m: update_status_callback(job['id'], m))
		except ZeroDivisionError:
			update_callback({"exception_note":"Bowtie returned no reads mapped to assembly; insert size not calculated.", 'lib_insert_length':'Unknown'})
	except subprocess.CalledProcessError as e:
		print job['job_type'], ": ", type(e), e, e.output
  		query("UPDATE assemblies SET status = 'exception', exception_note='{}:{}({})' WHERE id = '{}';".format(str(type(e)).encode('string_escape'), str(e).encode('string_escape'), str(e.output).encode('string_escape'), job['id']))
 	except ValueError as e:
 		#assembly exception
		print job['job_type'], ": ", type(e), e
  		query("UPDATE assemblies SET status = 'exception', exception_note='{}:{}' WHERE id = '{}';".format(str(type(e)).encode('string_escape'), str(e).encode('string_escape'), job['id']))
 	except KeyboardInterrupt:
 		query("UPDATE assemblies SET status='ready', exception_note='canceled by keyboard interrupt' WHERE id = '{}'".format(job['id']))
 		print "Terminated by keyboard interrupt."
 		quit()
 	except Exception as e:
 		print job['job_type'], ": ", type(e), e
 		if debug:
 			import traceback
 			traceback.print_exc(sys.stdout)
 			quit()
  		query("UPDATE assemblies SET status = 'exception', exception_note='{}:{}' WHERE id = '{}';".format(str(type(e)).encode('string_escape'), str(e).encode('string_escape'), job['id']))

	else:
		update_status_callback(job['id'], "Finished.") 
		query("""
UPDATE assemblies SET status='finished', 
					  average_coverage='{average_coverage}', 
					  num_contigs='{num_contigs}', 
					  n50='{n50}', 
					  num_bases='{num_bases}', 
					  assembly_version='{assembly_version}',
					  lib_insert_length='{lib_insert_length}',
					  dateCompleted='{date}',
					  fasta_file='{fasta_file}' 
					  WHERE id = '{id}';""".format(date=datetime.datetime.now().isoformat(),
												   id = job['id'],
												   **results))
	finally:
		if tempdir:
			shutil.rmtree(tempdir)
			
def main(loop=False):
	results = query("SELECT * FROM assemblies WHERE status='ready' LIMIT 1;")
	while True:
		while len(results) > 0:
			job = results[0]										 
			assemble(job, debug=('-debug' in sys.argv))
			results = query("SELECT * FROM assemblies WHERE status='ready' LIMIT 1;")
		
		if not loop:
			break
		sleep(30)
		
															 
if __name__ == "__main__":
	
	if '-test' in sys.argv:
		print "Starting test..."
		def query(s, f=lambda l: True):
			print "[SQL Cmd] {} ({})".format(s, f(0))
		job1 = {'id':'TEST', #fake job
				'job_type':'Supersembler',
				'path_1':'/shared/gn2/CFSANgenomes/CFSAN001656/CFSAN001656_01/',
				'file_1':'CFSAN001656_S8_L001_R1_001.fastq',
				'path_2':'/shared/gn2/CFSANgenomes/CFSAN001656/CFSAN001656_01/',
				'file_2':'CFSAN001656_S8_L001_R2_001.fastq',
				'ref_file':'',
				'ref_url':'',
				'trimmer':'null_trimmer',
				'trimmer_args':'',
				'insert_size':500,
				'min_contig_length':200,
				'k_value':177,
				'status':'ready',
				'exception_note':'',
				'dateAdded':datetime.datetime.today().strftime("%Y-%m-%d %h-%M-%S"),
				'dateCompleted':'',
				'accession':'CFSAN001656_01',
				'run_id':'',
				'fasta_file':'',
				'statter':'stat_fasta',
				'average_coverage':'20X',
				'num_contigs':'',
				'n50':'',
				'num_bases':'',
				'assembly_version':'',
				'lib_insert_length':''} 
		assemble(job1, debug=True)
		quit()
	
	
	#update options
	query("TRUNCATE TABLE assembly_options;")
	for (assembler_name, assembler_module) in assembler_dict.items():
		query("INSERT INTO assembly_options (option_type, option_value, option_description) VALUES ('assembler', '{}', '{} (Supports {})');".format(assembler_name, assembler_module.description.encode('string_escape'), ', '.join(assembler_module.supports)))
	for (trimmer, trim_func) in trimmer_dict.items():
		query("INSERT INTO assembly_options (option_type, option_value, option_description, option_supports) VALUES ('trimmer', '{}', '{}', 'All');".format(trimmer, str(trim_func.__doc__).encode('string_escape')))

	#run assemblies
	try:
		
		import prctl
		prctl.set_proctitle('Assem_dispatch')
		if '-nodaemon' in sys.argv:
			#raise ValueError
			main()
			quit()
		raise ValueError()
		import daemon
		print "Daemon module successfully imported; running in daemon mode."
		with daemon.DaemonContext(working_directory='/', detach_process=True, stdout=open('/data/Pipeline_jobs/assembly.log', 'a')):
			main(True)
	except (ImportError, ValueError) as e:
		print "Daemon module not found; running in normal mode."
		print e
		main(True)
	