import urllib
import MySQLdb
import json
import gzip

c = urllib.urlopen("http://smrt:8080/smrtportal/api/jobs/by-protocol/RS_HGAP_Assembly.1")

results = json.load(c)

c.close()

#parse results JSON objects and get ID's, then for each ID, retrieve assembly stats and load into DB

job_db = MySQLdb.connect(host="10.12.128.233", user="job_user", passwd="job_user", db="Jobs")
job = job_db.cursor(MySQLdb.cursors.DictCursor)

#
for i in [r['jobId'] for r in results['rows']]:
	m = json.load(urllib.urlopen("http://smrt:8080/smrtportal/api/jobs/{}/".format(i))) #m for message
	
	#print i, m.get("sampleName"), m.get('jobStatus'),
	
	if 'Completed' in m.get('jobStatus') and 'Multiple' not in m.get('sampleName'):
		job.execute("SELECT id FROM assemblies WHERE job_type LIKE '%{}';".format(i))
		if not job.fetchall: #uncollected assembly
			m2 = json.load(urllib.urlopen("http://smrt:8080/smrtportal/api/jobs/{}/metrics".format(i)))
			job.execute("""INSERT INTO assemblies (job_type,
												   status,
												   dateCompleted,
												   accession,
												   fasta_file,
												   average_coverage,
												   num_contigs,
												   n50,
												   num_bases,
												   assembly_version)
										   VALUES ('PacBio {
""".format(**m2)
			
	
	#print m.get('Job Name'), m.get('N50 Contig Length'), m.get('Sum of Contig Lengths'), m.get('Coverage')
	
	
	
job_db.close()