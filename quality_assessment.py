#!/usr/bin/env python

import os
import subprocess
import MySQLdb as mysqldb
import tempfile
import datetime
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import csv

"""

	Quality Assessment Daemon
	March 1, 2013
	Justin Payne
	
	MiSeq post-run quality assessment daemon, single-thread version. Interacts with Pipeline
	Job DB. Runs FASTX_Quality_Stats (FASTX Toolkit, Hannon Lab), then determines average
	read length and file size (in bytes).

"""


job_database = '/shared/gn3/job_tracking.db'
filestore = '/shared/gn2/CFSANgenomes/'

today = datetime.date.today().isoformat()


def qa_worker(job_id, fastq, job_db=None, callback=lambda s: None):

	if not job_db:
		job_db = mysqldb.connect(host='xserve15.fda.gov', user='job_user', passwd='job_user', db='Jobs')
	job = job_db.cursor()
	stats = fastq.replace('.fastq', '.stats').replace('.gz', '')
	plot = fastq.replace('.fastq', '.png').replace('.gz', '')
	summary = fastq.replace('.fastq', '.summary.stats').replace('.gz', '')
	q_plot = fastq.replace('.fastq', '.qscore.png').replace('.gz', '')
	
	job.execute("SELECT accession FROM qualities WHERE id = '{}';".format(job_id))
	job_accession = job.fetchone()[0]
	
	
	tempdir = ""
	
	try:
	
		if '.gz' in fastq:
			import gzip
			tempdir = tempfile.mkdtemp()
			fastq_basename = os.path.basename(fastq).replace('.gz', '')
			callback("Unzipping reads into {}".format(tempdir))
			with gzip.open(fastq) as r_in, open(os.path.join(tempdir, fastq_basename), 'w') as r_out:
				for l in r_in:
					r_out.write(l)
			fastq = os.path.join(tempdir, fastq_basename)
	
	
		with tempfile.NamedTemporaryFile() as temp_file:
			callback("QA started on {}".format(fastq))
			job.execute("UPDATE qualities SET status = 'running qa' WHERE id = '{}';".format(job_id))
			job_db.commit()
			callback("Quality-statting reads file")
			#print "fastx_quality_stats -Q 33 -i {fastq} -o {stats}".format(fastq=fastq, stats=stats)
			subprocess.check_call("fastx_quality_stats -Q 33 -i {fastq} -o {stats}".format(fastq=fastq, stats=stats), shell=True)
			callback("Running nucleotide distribution graph")
			subprocess.check_call("fastx_nucleotide_distribution_graph.sh -i {stats} -o {plot} -t {fastq}".format(fastq=fastq, stats=stats, plot=plot), shell=True)
			#print "awk '{{if(NR%4==2) print length($1)}}' {} > {}".format(fastq, temp_file.name)
			#print "awk 'BEGIN{{s=0;}}{{s=s+$1;}}END{{print s/NR;}}' {}".format(temp_file.name)
			callback("AWKing {}".format(temp_file.name))
			subprocess.check_call("awk '{{if(NR%4==2) print length($1)}}' {} > {}".format(fastq, temp_file.name), shell=True)
			result = subprocess.check_output("awk 'BEGIN{{s=0;}}{{s=s+$1;}}END{{print s/NR;}}' {}".format(temp_file.name), shell=True)
			#print result
			callback("Getting file size")
			file_size = os.path.getsize(fastq)
			num_reads = subprocess.check_output("grep -c '^@' {}".format(fastq), shell=True)
			try:
				import xmlrpclib
				s = xmlrpclib.ServerProxy('http://cfe1019692:8080')
				num_reads = int(num_reads)
				if num_reads < 37500:
					s.fire_event('QaqcFailureEvent', {'sample_id':job_accession, 'issues':['Insufficient yield ({} reads)'.format(num_reads),]})
					s.update_cfsan(s.get(job_accession)['Key'], 'QaIssues', s.get(job_accession)['QaIssues'] + "Insuff. reads ({})".format(num_reads))
			except:
				pass
		# Jamie's statting/plotting summary
		
		job.execute("UPDATE qualities SET status = 'running summary and plot' WHERE id = '{}';".format(job_id))
		job_db.commit()
		
		
		###Parse the output from fastx_quality_stats capturing column headers
		data = numpy.genfromtxt(stats, dtype='float',delimiter = '\t',skiprows=0, skip_header=0, skip_footer=0,usemask=False, names = True)

		####Calculate various summary statistics on the quality report
		cycles = max(data['column'])
		maxReads = max(data['count'])
		standardizedReadCounts = data['count']/maxReads
		avgReadLength = sum(standardizedReadCounts)
		avg_Q_score = sum(data['mean'])/cycles
		expectedCoverage = (sum(data['count']))/4500000

		callback("Writing stats summary")
		###Write quality summary stats to output file
		with open(summary, "w") as f:
			summary = csv.writer(f, delimiter = "\t",  quotechar="|", quoting=csv.QUOTE_MINIMAL)
			summary.writerow( ("CFSAN_ID ", "avgReadLength ", "avg_Q_score", "expectedCoverage") )
			summary.writerow( (os.path.splitext(stats)[0], avgReadLength, avg_Q_score, expectedCoverage) )
		f.close()


		callback("Plotting Q-score histogram")
		####Plot results
		matplotlib.pyplot.clf()
		matplotlib.pyplot.scatter((data['column']), (data['mean']))
		matplotlib.pyplot.ylim([0,40])
		matplotlib.pyplot.xlim([0,250])
		matplotlib.pyplot.ylabel('Mean Q-score')
		matplotlib.pyplot.xlabel('cycle number')
		matplotlib.pyplot.savefig(q_plot, format="png", dpi = 200) 

		job.execute("UPDATE qualities SET status = 'finished', average_read_length = '{}', date_completed = '{}',  file_size = {}, num_reads='{}' WHERE id = '{}';".format(str(result), today, str(file_size), str(num_reads), job_id))
		job_db.commit()
		callback("Finished.")
	except subprocess.CalledProcessError as e:
		callback("Subprocess error: {}".format(e))
		with open('/shared/gn3/job_logs/qa.log', 'a') as logfile:
			print e
			logfile.write('\nProblem with QA process on {} for job id {}. Error was:\n'.format(fastq, job_id))
			logfile.write(str(e))

			try:
				job.execute("UPDATE qualities SET status ='process_exception', exception='{}' WHERE id = '{}';".format(str(e).encode('string_escape'), job_id))
			except mysqldb.ProgrammingError:
				job.execute("UPDATE qualities SET status='unknown_exception', exception='Error could not be captured (probably due to formatting)' WHERE id = '{}';".format(job_id))
				job_db.commit()

	except Exception as e:
		callback("Error:{} {}".format(type(e), e))
		print stats, plot, summary, q_plot
		import traceback, sys
		traceback.print_exc(sys.stdout)
		try:
			job.execute("UPDATE qualities SET status='other_exception', exception='{}' WHERE id = '{}';".format(str(e).encode('string_escape'), job_id))
		except mysqldb.ProgrammingError as e:
			job.execute("UPDATE qualities SET status='unknown_exception', exception='Error could not be captured (probably due to formatting)' WHERE id = '{}';".format(job_id))
			print type(e), e
		job_db.commit()
		callback(e.filename)
		
	except KeyboardInterrupt as e:
		try:
			job.execute("UPDATE qualities SET status='ready' WHERE id = '{}';".format(job_id))
		except Exception:
			pass
		finally:
			raise e
	finally:
		try:
			os.remove(temp_file.name)
		except (OSError, UnboundLocalError):
			pass
		if tempdir:
			shutil.rmtree(tempdir)

if __name__ == "__main__":

	
	job_db = mysqldb.connect(host= 'xserve15.fda.gov', user='job_user', passwd='job_user', db='Jobs') #localhost
	job = job_db.cursor(mysqldb.cursors.DictCursor)

	job.execute("SELECT accession, file_path, file_root, id FROM qualities WHERE (status = 'ready' or status = 'priority') ORDER BY status;")
	
	for row in job.fetchall():
		fastq = os.path.join(filestore, row['file_path'], row['file_root'])
		def cb(s):
			print "[{}]\t{}".format(datetime.datetime.today().ctime(), s)
		#someday multiprocess this
		qa_worker(row['id'], fastq, job_db, cb)
		
	job_db.close()

