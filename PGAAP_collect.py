#!/usr/bin/env python

import MySQLdb as mysql
import ftplib
import os
import datetime
import shutil
import tarfile
import subprocess

pgaap_ftp = 'ftp-private.ncbi.nih.gov'


filestore = '/shared/gn2/CFSANGenomes/'
staging = '/shared/gn3/ncbi_submission/wgs_2013_pipeline/staging/'
returning = '/shared/gn3/ncbi_submission/wgs_2013_pipeline/returning/'
annotations = '/shared/gn3/FDAgenomes/'

usage = """

	PGAAP_collect daemon
	April 8, 2013
	Justin Payne
	
	Run as cron job, no arguments. Persistent data storage is a localhost MySQL db. Looks 
	for 'waiting' jobs in the db, checks to see if there's a tgz file waiting on the FTP 
	server, nabs it if so, unpacks it to the correct directory, and dispatches WGS
	submission job.


"""
try:
	print "Connecting to NCBI FTP server."
	ftp = ftplib.FTP(pgaap_ftp)
	ftp.login(user='fdacfsan', passwd='WmqwNb6L')
	print "Changing to AUTOSCAN directory."
	ftp.cwd(r"AUTOSCAN/")
except ftplib.Error as e:
	with open('/shared/gn3/ncbi_submission/wgs_2013_pipeline/wgs.log', 'a') as logfile:
		logfile.write("[{}]	FTP Error - couldn't connect to server. Error was:".format(datetime.datetime.today().isoformat()))
		logfile.write(str(e))
		raise e


job = mysql.connect(host='xserve15.fda.gov', user='job_user', passwd='job_user', db='Jobs') #xserve15

cur = job.cursor(mysql.cursors.DictCursor)




cur.execute("SELECT id, accession, group_tag, file_root FROM jobs WHERE status='waiting' AND job_type = 'PGAAP Submission'")


file_list = list()
ftp.dir(lambda d: file_list.append(d))
file_list = '\n'.join(file_list)

# print "Checking for completed PGAP runs:"
# for item in cur.fetchall():
# 	tarball = item['file_root'] + '.tgz'
# 	if tarball in file_list:
# 		print "Collecting {}...".format(tarball),
# 		with open('/shared/gn3/ncbi_submission/pgaap/pgaap_uploadReturn/{}.tgz'.format(item['file_root'], 'wb')) as dest_file:
# 			try:
# 				ftp.retrbinary('RETR {}'.format(tarball), dest_file.write)
# 				print "done."
# 				cur.execute("UPDATE jobs SET status = 'retrieved' WHERE accession = '{}' AND group_tag = '{}'  AND job_type = 'PGAAP Submission'".format(item['accession'], item['group_tag']))
# 				job.commit()
# 			except ftplib.Error as e:
# 				with open('/shared/gn3/ncbi_submission/wgs_2013_pipeline/wgs.log', 'a') as logfile:
# 					logfile.write('[{}]\tFTP Error: Failed to retreive {}. Error was'.format(datetime.datetime.today().isoformat(), tarball))
# 					logfile.write(str(e))
# 					raise e
# 			except IOError as e:
# 				with open('/shared/gn3/ncbi_submission/wgs_2013_pipeline/wgs.log', 'a') as logfile:
# 					logfile.write("[{}]	Filesystem Error: Couldn't save to {}. Error was".format(datetime.datetime.today().isoformat(), dest_file.name))
# 					logfile.write(str(e))
# 					raise e

print "Checking completed PGAP runs."
def process(arg, dirname, names):
	for node in names:
		if not "sqn" in node:
			return
		accession = node.split(".")[0]
		cur.execute("SELECT * FROM jobs WHERE accession='{}' AND job_type='PGAAP Submission' AND status NOT LIKE 'finished';".format(accession))
		results = cur.fetchall()
		errors = []
		if len(results):
			result = results[0]
			print "Collected an annotation for", accession
			#perform the action
			date_completed = datetime.date.today().isoformat()
			id = result['id']
			
			try:
				shutil.copyfile(os.path.join(dirname, node), os.path.join(filestore, result['file_path'] + ".sqn"))
				try:
					subprocess.check_call("tbl2asn -i {} -x .sqn -V b -o {}".format(
						os.path.join(os.path.join(dirname, node)
						os.path.join(filestore, result['file_path'] + ".gbk"), shell=True)
				except subprocess.CalledProcessError as e:
					errors.append(e)
				else:
					shutil.copyfile(os.path.join(filestore, result['file_path'] + ".gbk"),
									annotations)
				cur.execute("""UPDATE jobs SET date_completed='{date_completed}',
											   status='retrieved'
							   WHERE id='{id}';""".format(**locals()))
			except IOError as e:
				print "Couldn't remove file."
				print type(e), e
				errors.append(e)
			else:
				if not errors:
					os.remove(os.path.join(dirname, node))
			finally:
				job.commit()
			#test for completed group
			if result['group_tag']:
				cur.execute("SELECT DISTINCT status FROM jobs WHERE group_tag='{group_tag}' and group_tag NOT LIKE '';".format(**result))
				if len(cur.fetchall()) == 1:
					cur.execute("UPDATE jobs SET status='finished' WHERE group_tag='{group_tag}';".format(**result))
					print "Group {group_tag} finished! Closing out entries in jobs db.".format(**result)
					job.commit()
	if not len(os.listdir(dirname)):
		try:
			os.rmdir(dirname)
		except Exception as e:
			print "Couldn't clean up dir."
			print type(e), e
	
	
	
os.path.walk(returning, process, None)



# with open('/shared/gn3/ncbi_submission/pgaap/current_status.txt', 'w') as status_file:
# 	status_file.write('Current status of outstanding PGAAP jobs as of {}:'.format(datetime.date.today().strftime('%A, %d %B %Y %I:%M%p')))
# 	cur.execute("SELECT DISTINCT group_tag FROM jobs WHERE status NOT LIKE 'finished' AND job_type = 'PGAAP Submission';")
# 	for item in cur.fetchall():
# 		group_tag = item['group_tag']
# 		cur.execute("SELECT status FROM jobs WHERE group_tag = '{}' AND next_step = 'WGS Submission'".format(group_tag))
# 		rows = cur.fetchall()
# 		if all(map(lambda i: 'retrieved' in i['status'], rows)): #if all have been retrieved
# 			with open('/shared/gn3/ncbi_submission/wgs_2013_pipeline/wgs.log', 'a') as logfile:
# 				logfile.write('All annotations for {} have been downloaded and extracted as of {}'.format(group_tag, datetime.date.today().isoformat()))
# 				status_file.write('All annotations for {} complete.'.format(group_tag))
# 			#whole group has been retreived, so close it out as a job and run the subsequent steps
# 			cur.execute("UPDATE jobs SET status = 'finished', date_completed = '{}' WHERE group_tag = '{}'".format(datetime.date.today().isoformat(), group_tag))
# 			#move files to CFSAN filestore
# 			for row in cur.execute("SELECT accession, file_root FROM jobs WHERE group_tag = '{}';").fetchall():
# 				annot_dir = os.path.join(filestore, row['accession'], 'annot')
# 				if not os.path.exists(annot_dir):
# 					os.mkdir(annot_dir)
# 				for suffix in ('.fasta', '.email', '.template'):
# 					file = os.path.join(staging, row['file_root'] + suffix)
# 					shutil.move(file, '/shared/gn3/ncbi_submission/pgaap/pgaap_processed/')
# 				print "Extracting {} ...".format(row['file_root'] + '.tgz')
# 				tarball = tarfile.open(os.path.join(staging, row['file_root'] + '.tgz'), 'r:gz')
# 				tarball.extractall(annot_dir)
# 				
# 				cur.execute("INSERT INTO jobs (accession, file_root, job_type, status, date_added) VALUES ('{}', '{}', 'WGS Submission', 'ready', '{}');".format(row['accession'], row['file_root'], datetime.date.today().isoformat()))
# 		else:
# 			cur.execute("SELECT id FROM jobs WHERE group_tag = '{}' AND next_step = 'WGS Submission' AND status = 'retrieved'".format(group_tag))
# 			retrieved_rows = cur.fetchall()
# 			status_file.write('{} of {} finished for {}.'.format(len(retrieved_rows), len(rows), group_tag)
# 		job.commit()
# 	#	
		
		
#check for 'ready' submissions and upload
print "Uploading files:"
ftp.cwd("/")
cur.execute("SELECT id, accession, file_root FROM jobs WHERE status = 'ready' AND job_type = 'PGAAP Submission'")
for item in cur.fetchall():
	#upload waiting files
	try:
		with open(os.path.join(staging, item['file_root']), 'rb') as src_file:
			try:
				print "Uploading {}...".format(item['file_root']),
				ftp.storbinary('STOR {}'.format(item['file_root']), src_file)
				print "done."
				cur.execute("UPDATE jobs SET status = 'waiting' WHERE id = '{}'".format(item['id']))
			except ftplib.Error as e:
				with open('/shared/gn3/ncbi_submission/wgs_2013_pipeline/wgs.log', 'a') as logfile:
					logfile.write('[{}]\tFTP Error: Failed to upload {}. Error was'.format(datetime.datetime.today().isoformat(), src_file.name))
					cur.execute("UPDATE jobs SET status = 'upload_exception' WHERE id = '{}'".format(item['id']))
					logfile.write(str(e))
	except IOError as e:
		print "Didn't find {file_root}. Ignoring...".format(**item)


		
	job.commit()