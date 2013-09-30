#!/usr/bin/env python

import os
import subprocess
import MySQLdb as mysqldb
import tempfile
import datetime
import quality_assessment as qa

import multiprocessing


filestore = '/shared/gn2/CFSANgenomes/'

today = datetime.date.today().isoformat()



if __name__ == "__main__":
	job_db = mysqldb.connect(host='xserve15.fda.gov', user='job_user', passwd='job_user', db='Jobs') #localhost
	job = job_db.cursor(mysqldb.cursors.DictCursor)
	job.execute("SELECT accession, file_path, file_root, id FROM qualities WHERE (status = 'ready' or status LIKE '%_exception');")
	
	args = list()
	
	for row in job.fetchall():
		fastq = os.path.join(filestore, row['file_path'], row['file_root'])
		#load up args for pool processing
		args.append([row['id'], fastq])
	
	#for a in args:
		#print a
		
	#print len(args), " jobs."
		
	def worker(argv):
		qa.qa_worker(argv[0], argv[1])

	pool = multiprocessing.Pool(processes=multiprocessing.cpu_count() - 1) #8 cores, hyperthreaded
	pool.map(worker, args)
	
	job_db.close()

