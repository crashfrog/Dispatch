import MySQLdb as mysql
import assembly_dispatch
import multiprocessing
import datetime
import os

default_root = '/shared'

raise ValueError("Doesn't work. Some assemblers not thread-safe.")

def work(id):
	job_db = mysql.connect(host="xserve15.fda.gov", user='job_user', passwd='job_user', db='Jobs')
	job = job_db.cursor(mysql.cursors.DictCursor)
	j = job.execute("SELECT * FROM assemblies WHERE id='{}';".format(id))
	assembly_dispatch.assemble(j.fetchone(), job_db)

if __name__ == "__main__":
	job_db = mysql.connect(host="xserve15.fda.gov", user='job_user', passwd='job_user', db='Jobs')
	job = job_db.cursor(mysql.cursors.DictCursor)
	
	#update options
	job.execute("TRUNCATE TABLE assembly_options;")
	for (assembler_name, assembler_module) in assembly_dispatch.assembler_dict.items():
		job.execute("INSERT INTO assembly_options (option_type, option_value, option_description) VALUES ('assembler', '{}', '{}');".format(assembler_name, assembler_module.description.encode('string_escape')))
	for (trimmer, trim_func) in assembly_dispatch.trimmer_dict.items():
		job.execute("INSERT INTO assembly_options (option_type, option_value, option_description) VALUES ('trimmer', '{}', '{}');".format(trimmer, str(trim_func.__doc__).encode('string_escape')))
	job_db.commit()
	
	
	while (job.execute("SELECT COUNT(id), status FROM assemblies WHERE status='ready' GROUP BY status;")): #as long as there are some to get, get 8 or less and run on them
		core_load = 0
		ids = list()
		job.execute("SELECT id, job_type FROM assemblies WHERE status='ready' ORDER BY id LIMIT {};".format(multiprocessing.cpu_count()))
		for entry in job.fetchall():
			if core_load + assembly_dispatch.assembler_dict[entry['job_type']].core_load <= multiprocessing.cpu_count():
				job.execute("UPDATE assemblies SET status='claimed' WHERE id='{}';".format(entry['id']))
				ids.append(entry['id'])
				core_load += assembly_dispatch.assembler_dict[entry['job_type']].core_load
		job_db.commit()
		if ids:
			try:
				pool = multiprocessing.Pool(processes=8)
				pool.map(work, ids)
			except Exception:
				job.execute("UPDATE assemblies SET status='ready' WHERE id in ({});".format(",".join(["'{}'".format(i) for i in ids])))
				job_db.commit()