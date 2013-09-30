#!/usr/bin/env python

from Bio import SeqIO
import tempfile
import subprocess
import shutil
import os
import re


"""

	FASTA statting tool. Usable from the shell or imported as a module.

"""

entrez_CFSAN_genera = '("Campylobacter"[Organism]) OR ("Erwinia"[Organism]) OR ("Listeria"[Organism]) OR ("Escherichia"[Organism]) OR ("Vibrio"[Organism]) OR ("Salmonella"[Organism]) OR ("Bacillus"[Organism]) OR ("Achromobacter"[Organism]) OR ("Citrobacter"[Organism]) OR ("Proteus"[Organism]) OR ("Serratia"[Organism]) OR ("Brenneria"[Organism]) OR ("Paenibacillus"[Organism]) OR ("Brucella"[Organism]) OR ("Enterobacter"[Organism]) OR ("Clostridium"[Organism]) OR ("Cronobacter"[Organism]) OR ("Mycoplasma"[Organism]) OR ("Lymphocryptovirus"[Organism]) OR ("Klebsiella"[Organism]) OR ("Shigella"[Organism])'


def stat_fasta(fasta_file):
	"Basic FASTA statistics. Can't determine average coverage on its own."
	
	d = {'num_contigs':'',
		 'n50':'',
		 'num_bases':'',
		 'fasta_file':os.path.basename(fasta_file)
		 }
		 
	#print "statting ", fasta_file
	
	with open(fasta_file, 'r') as f:
		fasta = list(SeqIO.parse(f, "fasta"))
	
	if not len(fasta):
		raise ValueError("No contigs or improperly formatted fasta.")
		print "Tried to read:\n{}".format(open(fasta_file, 'rU').read())
		
	#determine number of contigs
	d['num_contigs'] = len(fasta)
	
	#determine N50, see http://en.wikipedia.org/wiki/N50_statistic
	f_prime = list()
	for contig in fasta:
		for i in range(len(contig)):
			f_prime.append(len(contig)) #"Create another list L' , which is identical to L, except that every element n in L has been replaced with n copies of itself"
	f_prime.sort()
	if len(f_prime) % 2 == 0: #is even:
		d['n50'] = (f_prime[len(f_prime) / 2] + f_prime[len(f_prime) / 2 + 1]) / 2
	else:
		d['n50'] = f_prime[len(f_prime) / 2]
		
	d['num_bases'] = sum([len(contig) for contig in fasta])
	
	return d
	
def stat_abyss(fasta_file, k_value=64):
	"FASTA statistics using information in ABySS headers, including average coverage."
	d = stat_fasta(fasta_file)
	
	contig_cov = list()
	
	with open(fasta_file, 'r') as f:
		for contig in SeqIO.parse(f, "fasta"):
			raw_cov = float(contig.description.split(" ")[2])
			contig_cov.append(raw_cov / (float(len(contig)) - float(k_value) + 1))
	
	avg_cov = sum(contig_cov) / float(len(contig_cov))
	d['average_coverage'] = "{}X".format(int(avg_cov))
	
	return d
	
def stat_velvet(fasta_file, k_value=171):
	"FASTA statistics using information in Velvet headers, including average coverage."
	d = stat_fasta(fasta_file)
	
	#print "statting (velvet headers)", fasta_file
	
	#>NODE_1_length_71394_cov_21.306412
	cov_exp = re.compile(r"length_(?P<length>\d*)_cov_(?P<cov>\d*\.\d*)")
	
	with open(fasta_file, 'r') as f:
		covs = list()
		for contig in SeqIO.parse(f, "fasta"):
			m = cov_exp.search(contig.description)
			if m:
				cov=(float(m.groupdict()['cov']))
				length=(float(m.groupdict()['length']))
				covs.append(cov * (length - int(k_value) + 1) / length)
						
	avg_cov = sum(covs) / float(len(covs))
	d['average_coverage'] = "{}X".format(int(avg_cov))
	
	return d
		
def stat_blast(fasta_file, callback=None, organism=entrez_CFSAN_genera, num_threads=4):
	"Experimental method to determine assembly 'realness' by blasting against reference sequences."
	print ""

	import datetime

	if not callback:
		def callback(s):
			print "[{}]".format(datetime.datetime.today().ctime()), s

	callback("Importing modules...")

	from Bio.Blast import NCBIWWW
	import xml.etree.ElementTree as xml
	
	callback("Statting fasta...")
	d = stat_fasta(fasta_file)
	
	results = list()
	
	callback("Loading {}...".format(fasta_file))
	
# 	with open(fasta_file, 'r') as f:
# 		os.chdir("/data/blast_db")
# 		for contig in sorted(list(SeqIO.parse(f, "fasta")), key=lambda c: -len(c)):
# 			if len(contig) >= d['n50']:
# 				callback("BLASTing {} ({} bases)...".format(contig.description, len(contig)))
# 				result = xml.parse(subprocess.check_output("blastn -db refseq_genomic -query {} -task blastn -num_threads 8 -outfmt 5 -max_target_seqs 1".format(contig.seq), shell=True)
# 				#results.append(result.findall(".//Hsp_identity")[0:1])
# 				realness = float(result.findall(".//Hsp_identity")[0].text) / float(len(contig)) * 100.0
# 				print "{}% identity to something real".format(realness)

				
	with open(fasta_file, 'r') as f:
	
		os.chdir("/data/blast_db")
	
		#contigs = filter(lambda c: len(c) >= d['n50'], sorted(list(SeqIO.parse(f, "fasta")), key=lambda c: -len(c)))
		contigs = filter(lambda c: len(c) < 10000, sorted(list(SeqIO.parse(f, "fasta")), key=lambda c: -len(c)))
		p = 0
		results = list()
		while p < len(contigs):
			q = p + num_threads
			
			
			#make a query file
			
			with tempfile.NamedTemporaryFile("w") as query_file:
				query = "\n".join([">{}\n{}".format(c.description, c.seq) for c in contigs[p:q]])
				query_file.write(query)
				callback("Made contigs file {}".format(query_file.name))
				callback("BLASTing {}-{} of {} contigs...".format(p + 1, min(q, len(contigs)), len(contigs)))
				
				r = subprocess.check_output("blastn -db refseq_genomic -query {} -task blastn -num_threads 8 -outfmt 5 -max_target_seqs 1".format(query_file.name), shell=True)
				callback("BLAST complete. Parsing...")
				result = xml.parse(r)
				#result.write("/home/justin.payne/blast{}.xml".format(min(q, len(contigs))))
				results.append(result)
				p = q
		callback("Done.")
		realnesses = list()
		for r in results:
			for iteration in r.findall(".//Iteration"):
				realness = float(iteration.find("/Hit/Hit_hsps/Hsp/Hsp_identity").text) / float(iteration.find("/Iteration_query-len").text) * 100.0
				realnesses.append(realness)
				
		total_realness = sum(realnesses) / float(len(realnesses))
		
		return total_realness
		
def find_closest_ref(fasta_file, callback=None, update_callback=lambda d: None, organism=entrez_CFSAN_genera):
	"Find closest match in NCBI Refseq to longest contig, then collect URL for it"
	if not callback:
		import datetime
		def callback(s):
			print "[{}]".format(datetime.datetime.today().ctime()), s
			
	callback("Importing modules...")

	from Bio.Blast import NCBIWWW
	import xml.etree.ElementTree as xml
	
	callback("Loading fasta ({})...".format(fasta_file))
	with open(fasta_file, 'r') as f:
		contigs = iter(sorted(list(SeqIO.parse(f, 'fasta')), lambda a,b: cmp(len(a), len(b))))
		contig = contigs.next()
		while len(contig) < 1500:
			try:
				contig = contigs.next()
			except StopIteration:
				break
	callback("Longest contig is {} bases. BLASTing...".format(len(contig)))
	r = NCBIWWW.qblast("blastn", "chromosome", ">{}\n{}".format(contig.description, contig.seq), 
					   alignments=1, 
					   entrez_query="{}".format(organism),
					   hitlist_size=1,
					   filter='L')
	callback("BLAST finished.")
	result = xml.parse(r)
	refseq = result.find(".//Iteration/Iteration_hits/Hit/Hit_id").text.split("|")[1]
	refseq_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text'.format(refseq)
	update_callback({'ref_file':refseq, 'ref_url':refseq_url})
	return refseq



# 				
# 	with open(fasta_file, 'r') as f:
# 		contigs = filter(lambda c: len(c) >= d['n50'], sorted(list(SeqIO.parse(f, "fasta")), key=lambda c: -len(c)))
# 		p = 0
# 		results = list()
# 		while p < len(contigs):
# 			q = p + num_threads
# 			callback("BLASTing {}-{} of {} contigs...".format(p, min(q, len(contigs)), len(contigs)))
# 			r = NCBIWWW.qblast("blastn", "refseq_genomic", "\n".join([">{}\n{}".format(c.description, c.seq) for c in contigs[p:q]]), alignments=1, entrez_query="{}[Organism]".format(organism))
# 			callback("BLAST complete. Parsing...")
# 			result = xml.parse(r)
# 			result.write("/home/justin.payne/blast{}.xml".format(min(q, len(contigs))))
# 			results.append(result)
# 			p = q
# 		callback("Done.")
		
	
def quast_compare(path, fastas, gi=None, callback=None, update_callback=lambda d: None, debug=False, **kwargs):
	"Find best reference, then use Quast to compare assemblies and return best one."
	if not callback:
		import datetime
		def callback(s):
			print "\n[{}] ".format(datetime.datetime.today().ctime()) + s
			
	
	import urllib
	import tempfile
	import subprocess
	import glob
	import csv
	def rank(a, b):
		"Rank two assemblers based on a priority list"
		#our own intuition about which assemblies are better, from best to worst
		assemblers = ['SPAdes', 'ABySS', 'Velvet', 'CLC']
		if a['assembler'] in assemblers and b['assembler'] in assemblers:
			return cmp(assemblers.index(a['assembler']),
					   assemblers.index(b['assembler']))
		elif a['assembler'] in assemblers:
			return -1
		else:
			return 1
	temp_dir = tempfile.mkdtemp()	
	def cb(a, b, c):
		if int(a) % 100 == 0:
			callback("Downloaded block {} of {}\r".format(a, int(c) / int(b)))
	try:	
		fastas.sort(rank)
		if "WORST" in fastas[0]['assembler']:
			raise ValueError("MISASSEMBLY DO NOT USE")
		callback("Getting closest reference.")
		if not gi:
			gi = find_closest_ref(os.path.join(path, fastas[0]['fasta_file']), callback=callback, update_callback=update_callback, **kwargs)
		urllib.urlretrieve('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text'.format(gi),
						   os.path.join(temp_dir, "reference.fasta"),
						   cb)
		callback("Getting gene annotations.")
		urllib.urlretrieve('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=ft&retmode=text'.format(gi),
						   os.path.join(temp_dir, "reference.genes"),
						   cb)
		with open(os.path.join(temp_dir, "reference.fasta"), 'r') as ref_fasta:
			update_callback({'ref_url':ref_fasta.readline()})
			
		if os.path.exists("{}quast/".format(path)):
			shutil.rmtree("{}quast/".format(path))
			
		for r in fastas:
			if not os.path.exists(os.path.join(path, r['fasta_file'])):
				callback("Symlinking {}".format(r['fasta_file']))
				try:
					os.symlink(os.path.join(path, r['fasta_file'].split('.')[0]+".fasta"), os.path.join(path, r['fasta_file']))
				except OSError:
					try:
						shutil.copyfile(os.path.join(path, r['fasta_file'].split('.')[0]+".fasta"), os.path.join(path, r['fasta_file']))
					except OSError:
						pass
		callback("Running QUAST with reference gi{}...".format(gi))
		subprocess.check_call("quast {} -R {}/reference.fasta -G {}/reference.genes -o {} --labels {} --gene-finding -t 8".format(
							  " ".join([os.path.join(path, r['fasta_file']) for r in fastas]),
							  temp_dir,
							  temp_dir,
							  os.path.join(path, "quast"),
							  ",".join([r['assembler'] for r in fastas])), shell=True)
	except OSError:
		pass
	finally:
		callback("Cleaning {}...".format(temp_dir))
		shutil.rmtree(temp_dir)

	
	with open("{}/quast/report.tsv".format(path), 'rU') as report:
		
		#load generated report
		r = list(csv.DictReader(report, dialect='excel', delimiter='\t'))
		#QUAST reports come in with the assemblers as the columns; need to pivot this table so assemblers are as rows
		tbl = [] #table
		#get assemblers from column headers
		for h in r[0].keys():
			if 'Assembly' not in h:
				tbl.append({'assembler':h})
		for row in r:
			#these are the specific fields we want, in a row-dict whose key is "Assembly" (I know, right?)
			for h in row.keys():
				if 'Assembly' not in h:
					tr = filter(lambda t: h in t['assembler'], tbl)[0]
					tr[row['Assembly']] = row[tr['assembler']]
		for r in tbl:
			for h in r.keys():
				#cast to number values if possible
				try:
					r[h] = int(r[h])
				except ValueError:
					try:
						r[h] = float(r[h])
					except ValueError:
						if r[h] == '-':
							r[h] = 0
			for h in ('Genome fraction (%)', 'NA50', 'Reference length', 'Reference GC (%)', '# misassemblies'):
				if h not in r:
					r[h] = 0
			try:
				#"goodness" model
				r['goodness'] = ((r['Genome fraction (%)'] / 100) or 0) * ( (r['NA50'] or r['N50']) * min(r['Total length'], r['Reference length']) )  / ( pow( r['GC (%)'] - r['Reference GC (%)'], 2) * max(r['# misassemblies'], 1) * r['# contigs'])
			except KeyError as e:
				if debug:
					print tbl
					raise e
		
	return tbl
		
	
	
	
	
	
def read_map(assembly, read1, read2=None, callback=lambda s: None):
	"Remap raw reads to assembly using Bowtie to determine average library insert length."

	temp = tempfile.mkdtemp()
	index = os.path.join(temp, "index")
	#alignment = os.path.join(temp, "alignment.map")
	try:
		callback("statting: building index")
		results = subprocess.check_output("bowtie-build {assembly} {index}".format(assembly=assembly, index=index), shell=True)
		callback("statting: mapping reads to assembly")
		if read2:
			results = subprocess.check_output("bowtie -p 4 -v 1 -s 25 {index} -1 {read1} -2 {read2} -I 50 -X 3000 --quiet --suppress 1,2,6,7,8".format(index=index, read1=read1, read2=read2), shell=True)
		else:
			results = subprocess.check_output("bowtie -p 4 -v 1 -s 25 {index} -1 {read1} --quiet --suppress 1,2,6,7,8".format(index=index, read1=read1), shell=True)
		
		#capture results
		
		results = results.split("\n") #turn into lines
		inserts = list()
		for i in range(0, len(results), 2):
			try:
				if results[i].split("\t")[0] not in results[i+1].split("\t")[0]:
					print results[i].split("\t")[0], results[i+1].split("\t")[0]
					raise ValueError("Paired-end reads didn't map to same contigs.")
				coord1 = float(results[i].split("\t")[1])
				coord2 = float(results[i+1].split("\t")[1]) + float(len(results[i+1].split("\t")[2]))
				inserts.append(coord2 - coord1)
			except IndexError:
				pass
				
		avg_insert = sum(inserts) / float(len(inserts))
		

		
	finally:
		shutil.rmtree(temp)
		
	return avg_insert
		
#fasta_statter.read_map("/shared/gn2/CFSANgenomes/CFSAN002091/asm/CFSAN002091_01.fasta", "/shared/gn2/CFSANgenomes/CFSAN002091/CFSAN002091_01/CFSAN002091-01_S6_L001_R1_001.fastq", "/shared/gn2/CFSANgenomes/CFSAN002091/CFSAN002091_01/CFSAN002091-01_S6_L001_R2_001.fastq")
#fasta_statter.stat_blast("/shared/gn2/CFSANgenomes/CFSAN002091/asm/CFSAN002091_01.fasta")
	
	
if __name__ == "__main__":
	import sys
	import csv
	import traceback
	def ucb(d):
			for (k,v) in d.items():
				print k, ':', v
	debug = '-debug' in sys.argv
	if '-test-quast' in sys.argv:
		fastas = list(csv.DictReader(open('/shared/gn2/CFSANgenomes/CFSAN004357/asm/comparative_assembly_stats.txt', 'r'), delimiter='\t'))
		print quast_compare('/shared/gn2/CFSANgenomes/CFSAN004357/asm/',
					  fastas,
					  #gi = '523917454',
					  update_callback=ucb,
					  debug=debug)
	elif '-redo-quast' in sys.argv:
		try:
			accession = sys.argv.pop(sys.argv.index('-redo-quast') + 1)
		
			fastas = list(csv.DictReader(open('/shared/gn2/CFSANgenomes/{}/asm/comparative_assembly_stats.txt'.format(accession), 'r'), delimiter='\t'))
			quast_results = quast_compare('/shared/gn2/CFSANgenomes/{}/asm/'.format(accession),
						  fastas,
						  update_callback=ucb,
						  debug=debug)
			headers = ('assembler',
				   'n50',
				   'num_contigs',
				   'Total length',
				   'Reference length',
				   'Genome fraction (%)',
				   '# misassemblies',
					'goodness')
		
			widths = [max([max(len(str(r.get(h, ''))), len(str(h))) for r in quast_results]) + 2 for h in headers] #get longest string in each field in quast_results
			print "".join([str(h).ljust(w) for (w, h) in zip(widths, headers)])
			for r in sorted(quast_results, key=lambda r: r['assembler']):		
				print "".join([str(r.get(h, '-')).ljust(w) for (w, h) in zip(widths, headers)])
		except IOError as e:
			print "{}: Supersembler not performed on this isolate.".format(e)
			try:
				subprocess.check_call("quast /shared/gn2/CFSANgenomes/{0}/asm/*.fasta -o /shared/gn2/CFSANgenomes/{0}/asm/quast/ -R /shared/gn2/CFSANgenomes/CFSAN002060/asm/CFSAN002060.fasta --gene-finding -t 8".format(accession), shell=True)
			except subprocess.CalledProcessError as e:
				print e.output
		except IndexError:
			print "Specify FDA CFSAN accession number (i.e. 'CFSAN001250')"
			quit()
	elif '-scan-assemblies' in sys.argv:
		import glob
		for path in glob.glob('/shared/gn2/CFSANgenomes/*/asm/'):
			try:
				fastas = list(csv.DictReader(open(os.path.join(path, 'comparative_assembly_stats.txt'), 'r'), delimiter='\t'))
				quast_compare(path, fastas, update_callback=ucb)
			except IOError:
				try:
					subprocess.check_call("quast {0}*.fasta -o {0}quast/ --gene-finding -t 8".format(path), shell=True)
				except subprocess.CalledProcessError:
					pass
			except Exception:
				traceback.print_exc(file=sys.stderr)
	else:	
		ucb(stat_fasta(os.path.normpath(sys.argv[1])))
	