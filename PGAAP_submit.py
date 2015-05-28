#!/usr/bin/env python

import ODBCDictionary
import CFSANPreferences
import pyodbc
import datetime
import MySQLdb as mysqldb
import shutil
import glob
import os
import re

usage = """
PGAAP_submit April 8, 2013 by Justin Payne

Prepare isolate assemblies for submission to PGAAP, copy to staging area, and upload to the PGAAP server,
then update the Genomics database and PGAAP_collect job database. First argument is an arbitrary
group tag, this is checked for uniqueness in the DB and used to group individual submissions together
to track completeness. No spaces, but it doesn't matter what you put in here.

Usage examples:
	
	PGAAP_submit april_cbots CFSAN001992 CFSAN001993 CFSAN001994...
	
	or
	
	PGAAP_submit bobs_heidelbergs submission_metadata.txt
	
	Note - you can use a tab or comma-delimited metadata sheet, but PGAAP_submit ignores
	everything but the first column (where it expects to find CFSAN accession numbers.)
	
	This module can also be imported into BioNumerics if you have the CFSAN BioNumerics
	python modules.
	
"""

emails = """Errol.Strain@fda.hhs.gov
Justin.Payne@fda.hhs.gov
Yan.Luo@fda.hhs.gov
Ruth.Timme@fda.hhs.gov"""

template_path = 'X:/ncbi_submission/pgaap/submit.template' #change these to preference items, later

filestore = CFSANPreferences.prefs['CFSAN genomes directory']

staging = 'X:/ncbi_submission/pgaap/staging/'


def submit(isolate_ids, group_tag, callback=None, continue_to_wgs=True, copy_fasta=False):

	next_step = ''
	if continue_to_wgs:
		next_step = 'WGS Submission'
		
	db = pyodbc.connect(CFSANPreferences.prefs["Genomics ODBC connection string"])
	#job = sqlite3.connect(job_database, autocommit=False)
	job = mysqldb.connect('cfu5059592', 'django', 'django', 'Jobs')

	isolate_list = ODBCDictionary.ODBCFactoryById(isolate_ids, db=db)
	currDate = datetime.date.today().strftime('%Y%m%d')
	template = open(template_path, 'rU').read()
	
	#basic sanity checks - is the required information in the DB? Is there actually an assembly file to submit?
	
	#check for presence of required fields:
	missing_list = list()
	for v in ('NCBI_LocusTag', 'NCBI_TaxID', 'NCBI_BioProject', 'StrainName', 'AssemblyVersion', 'PeakDepth_AvgCoverage'):	
		for i in isolate_list:
			if not i[v]:
				missing_list.append('{0} missing {1}.'.format(i['Key'], v))
	if len(missing_list) > 0:
		raise ValueError('Missing information:\n' + '\n'.join(missing_list))
																												   
	#check that there's actually an assembly file:	
	for i in isolate_list:
		if not len(glob.glob(os.path.join(filestore, i.FdaAccession, 'asm', i.FdaAccession + '.f*a'))) > 0:
			missing_list.append('No assembly fasta or fna present in {}.'.format(os.path.join(filestore, i.FdaAccession, 'asm')))
	if len(missing_list) > 0:
		raise ValueError('Missing assembly file or files:\n' + '\n'.join(missing_list))
	
	
	if not callback:
		def callback(x): #do-nothing
			pass
	
	for isolate in isolate_list:
		#
		file_root = '{}_{}'.format(isolate.NCBI_LocusTag, currDate)
		print "File root: ", file_root
		
		with open(os.path.join(filestore, isolate.FdaAccession, isolate.FdaAccession + '.meta'), 'w') as metafile:
			for (k, v) in CFSANPreferences.metafile_keys.iteritems():
				try:
					metafile.write('{0}={1}'.format(k, v(isolate)))
				except TypeError: #item was not a function
					metafile.write('{0}={1}'.format(k, isolate[v]))
					
		with open(os.path.join(filestore, isolate.FdaAccession, 'annot', file_root + '.cmt'), 'w') as structured_comment_file:
			for (k, v) in CFSANPreferences.structured_comment_keys.iteritems():
				try:
					structured_comment_file.write('{0}={1}'.format(k, v(isolate)))
				except TypeError: #item was not a function
					structured_comment_file.write('{0}={1}'.format(k, isolate[v]))
					
		#build email file and submission template
		with open(os.path.join(staging, file_root + '.email'), 'w') as email_file:
			print "Writing email file."
			email_file.write(emails)
			
		with open(os.path.join(staging, file_root + '.template'), 'w') as template_file:
			print "Writing template file."
			template_file.write(template)
		
		#build contig header
		header_line = ''
		for (k, v) in CFSANPreferences.pgaap_contig_keys.iteritems():
			try:
				value = v(isolate)
			except TypeError:
				value = isolate[v]
			if value:
				header_line = header_line + "[{0}={1}]".format(k, value)
		header_line = header_line + "[gcode=11]"
		print "Header: ", header_line
			
		#get fasta files
		if not copy_fasta:
			with open(glob.glob(os.path.join(filestore, isolate.FdaAccession, 'asm', isolate.FdaAccession + '.f*a'))[0], 'r') as src_fasta, open(os.path.join(staging, file_root + '.fasta'), 'wb') as dest_fasta:
				print "Writing fasta file."
				contig_number = 0
				for currLine in src_fasta:
					if '>' in currLine[0:1]:
						dest_fasta.write('>gn1|fda|{}_{:05} {}\n'.format(isolate.NCBI_LocusTag, contig_number, header_line))
						contig_number = contig_number + 1
					else:
						#do we need to fix it? exclude short contigs?
						dest_fasta.write(currLine)
		else:
			shutil.copy(glob.glob(os.path.join(filestore, isolate.FdaAccession, 'asm', isolate.FdaAccession + '.f*a'))[0], os.path.join(staging, file_root + '.fasta'))
		
		#edit headers
		
		#update db's
		print "Updating databases and making job."
		job.cursor().execute("INSERT INTO jobs (accession, group_tag, job_type, status, file_path, file_root, next_step, date_added) VALUES ('{}', '{}', 'PGAAP Submission', 'ready', '{}', '{}', '{}', '{}')".format(isolate.Key, group_tag, isolate.NCBI_LocusTag, file_root, next_step, datetime.date.today().isoformat()))
#		job.cursor().execute("""INSERT INTO metadata (accession,  locus_tag,  bioproject,  biosample,  taxid,  full_name,  subspecies,  serovar,  strain,  isolate,
#			    		specific_host,  county,  state,  isolation_date,  isolation_source,  culture_collection,  assembler_version,  depth,  sequencing_tech) 
#			    	   VALUES ('{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}')""".format(isolate.Key,
#			    	   																		 isolate.NCBI_LocusTag,
#			    	   																		 isolate.NCBI_BioProject,
#			    	   																		 isolate.NCBI_BioSample,
#			    	   																		 isolate.NCBI_TaxID,
#			    	   																		 ODBCDictionary.FullName(isolate),
#			    	   																		 isolate.Subspecies,
#			    	   																		 isolate.Serovar,
#			    	   																		 isolate.StrainName,
#			    	   																		 isolate.NCBI_SpecificHost,
#			    	   																		 isolate.Country,
#			    	   																		 isolate.State,
#			    	   																		 isolate.IsolationDate,
#			    	   																		 isolate.IsolationSource,
#			    	   																		 isolate.NCBI_CultureCollection,
#			    	   																		 isolate.AssemblyVersion,
#			    	   																		 isolate.PeakDepth,
#			    	   																		 isolate.ASMSequenceTechnology))
#		db.cursor().execute("UPDATE GenomicsBN.dbo.ENTRYTABLE SET Annotation = ? WHERE [KEY] = ?", ('submitted to PGAAP {}'.format(datetime.date.today().isoformat()), isolate.Key))
		job.commit()
		db.commit()
		#
	db.close()
	job.close()
	
	
	
if __name__ == '__main__':
	#invoked from command line with either file as sole argument or list of cfsan ids
	import sys
	try:
		#argument vector may be a filename in the cwd, try to open it
		metadata_file = open(sys.argv[2], 'rU')
		#sniff file and try to determine delimiter and presence of headers
		import csv
		sample = metadata_file.read(1024)
		metadata_file.seek(0) #return to beginning of file
		dialect = csv.Sniffer().sniff(sample)
		if csv.Sniffer().has_header(sample):
			metadata = csv.DictReader(metadata_file, dialect=dialect)
		else:
			metadata = csv.reader(metadata_file, dialect)
		ids = map(lambda l: l[0], metadata)
		
		
	except IOError:
		if 'CFSAN' in sys.argv[2]:
			ids = sys.argv[2:] #arg vector is list of cfsan ids
		else:
			print usage
			raise ValueError('Invalid argument - provide group tag and filename or space-delimited CFSAN ids')
	except IndexError: #invoked with no arguments
		print usage
		quit()
	finally:
		submit(ids, sys.argv[1])