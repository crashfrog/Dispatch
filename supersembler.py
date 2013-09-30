#!/usr/bin/env python
#import assembly_dispatch
import subprocess
import shutil
import os
import fasta_statter

"""
	Supersembler
	
	Meta-assembly binding for Pipeline: Assembly Dispatch. Calls all the assemblers in
	order, compares the results to find the assembly with the lowest contigs, and returns
	it. Future improvements: Use the CGAL assembly likelihood measure instead of simple
	contig number (via Rahman and Pachter 2013).

"""

description = "Justin's round-robin assembler shoot-out. Runs them all and selects the best assembly."

core_load = 8

dont_use = ("Supersembler", "Velvet_optimize")

supports = ('MiSeq', 'IonTorrent')

def assemble(assembler_dict, data_type, ref_file=None, callback=lambda s: None, update_callback=lambda s: None, debug=True, **kwargs):

	def update_callback_hook(d):
		"""Hook function to capture statting information."""
		pass 

	results = list()
	errors = list()
	
	kwargs['debug'] = debug
	
	for (assembler_name, assembler) in assembler_dict.items():
		if assembler_name not in dont_use: #infinite loop danger
			def callback_wrapper(s): #closure to wrap callback method, if any, from assembly_dispatch
				callback("{}: {}".format(assembler_name, s))
		
			parameters = {'data_type':data_type}
			parameters.update(kwargs) # copy parameter keys/values
			parameters['callback'] = callback_wrapper
			#del parameters['update_callback']
			parameters['fasta_file_name'] = "{}.{}.fasta".format(parameters['accession'], assembler_name)
			try:
				if data_type in assembler.supports:
					callback("Starting {}...".format(assembler_name))
					r = assembler.assemble(**parameters)
					r['assembler'] = assembler_name
					results.append(r)
			except subprocess.CalledProcessError as e:
				print assembler_name, ": ", type(e), e, e.output
				errors.append("{}:{}{} ({})\n".format(assembler_name, type(e), e, e.output))
			except Exception as e:
				import traceback
				import sys
				if debug:
					traceback.print_exc(file=sys.stderr)
				callback("{}:{}{}\n".format(assembler_name, type(e), e))
				errors.append("{}:{}{}\n".format(assembler_name, type(e), e))
	
	#determine which is best
	f = "{}.fasta".format(kwargs['accession'])
		
	callback("{accession}: running comparison".format(**kwargs))
	
	try:	
		d = results[0]
	except (IndexError, ValueError) as e:
		errors.append(e)
		raise ValueError("No results from assembly. Errors were:{}".format(errors))
	
	
	#print "completed assemblies: ", len(results)
	import csv
	callback("Writing {}...".format(os.path.join(kwargs['path'], "comparative_assembly_stats.txt")))
	with open(os.path.join(kwargs['path'], "comparative_assembly_stats.txt"), 'w') as statsfile:
		headers = sorted(results[0].keys())
		dw = csv.DictWriter(statsfile, fieldnames=headers, extrasaction='ignore', delimiter='\t')
		#print results[0].keys().sort()
		dw.writeheader()
		dw.writerows(results)
	#print a nice visual table
	widths = [max([max(len(str(r.get(h, ''))), len(str(h))) for r in results]) + 2 for h in headers] #get longest string in each field in results
	print "".join([str(h).ljust(w) for (w, h) in zip(widths, headers)])
	for r in results:		
		print "".join([str(r.get(h, '-')).ljust(w) for (w, h) in zip(widths, headers)])
	try:
	
		quast_results = fasta_statter.quast_compare(kwargs['path'], results, callback=callback, update_callback=update_callback, gi=ref_file, debug=debug)	
		
		# === determine best assembly ===
		#generate "goodness" score
		

		
		for r in quast_results:
			r.update(filter(lambda rs: rs['assembler'] in r['assembler'], results)[0]) #merge two results
			try:
				if r['goodness'] > d['goodness']:
					d = r
			except KeyError:
				if r['n50'] > d['n50']:
					d = r
		#headers = sorted(quast_results[0].keys())
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
			
			
		
		
			
			
	except Exception as e:
		callback("Quast module didn't complete, comparing N50. Error was {}:{}".format(str(type(e)), str(e)))
		if debug:
			import traceback
			class CallbackWrapper(object):
				"duck-typed object to route stack traceback through callback system"
				def write(s):
					for line in s.split("\n"):
						callback(line)
			traceback.print_exc(CallbackWrapper())
		for r in results:
			if (r['n50'] > d['n50']) and 'WORST' not in r['assembler']:
				d = r
	
	try:
		shutil.move(os.path.join(kwargs['path'], d['fasta_file']), os.path.join(kwargs['path'], f))
	except shutil.Error:
		try:
			shutil.move(os.path.join(kwargs['path'], "{}.{}.fasta".format(parameters['accession'], d['assembler'])), os.path.join(kwargs['path'], f))
		except shutil.Error:
			print "Couldn't find assembly " + d['fasta_file']
	d['fasta_file'] = f
	
	try:
		os.remove("{path}/purposeful_misassembly.fasta".format(**kwargs))
	except OSError:
		pass
	return d
	
if __name__ == "__main__":
	#debug
	import datetime
	import assembly_dispatch
	def cb(d):
		print "[{}] {}".format(datetime.datetime.today().ctime(), d)
	def bcb(d):
		for (k, v) in d.items():
			cb("{} : {}".format(k, v))
	d = {}
	d.update(assembly_dispatch.assembler_dict)
	del d['SPAdes'] #when spades dies it kills the Python interpreter
	del d['WORST'] #waste of time in debugging
	print assemble(assembler_dict=d,
			 path='/home/justin.payne/', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001656/CFSAN001656_01/CFSAN001656_S8_L001_R1_001.fastq', 
			 reads2='/shared/gn2/CFSANgenomes/CFSAN001656/CFSAN001656_01/CFSAN001656_S8_L001_R2_001.fastq', 
			 accession='CFSAN001656_01', 
			 insert_size=500,
			 callback=cb,
			 update_callback=bcb,
			 k_value=143,
			 ref_file='528816716',
			 debug=True)
