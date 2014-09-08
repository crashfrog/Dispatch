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
	job.execute("SELECT accession, file_path, file_root, id FROM qualities WHERE (status = 'ready' or status = 'priority') ORDER BY status;")
	
	args = list()
	
	for row in job.fetchall():
		fastq = os.path.join(filestore, row['file_path'], row['file_root'])
		#load up args for pool processing
		args.append([row['id'], fastq])
	
	#for a in args:
		#print a
		
	#print len(args), " jobs."
		
	def worker(argv):
		try:
			qa.qa_worker(argv[0], argv[1])
		except Exception as e:
			import traceback
			traceback.print_exc()

	pool = multiprocessing.Pool(processes=16) #8 cores, hyperthreaded, 30% load roughly
	pool.map(worker, args)
	
	job_db.close()

