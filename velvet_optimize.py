#!/usr/bin/env python

import subprocess
import shutil
import os
import tempfile
import datetime
import fasta_statter


import velvet_assembler_notrimmer
import worst_assembler

"""
	
	Velvet "Optimizer"-style determination of local-optimum kmer length for Assembly
	Dispatch.

"""

description = "Two-round local-maximum k-value optimization assembly using Velvet."

core_load = 1

supports = velvet_assembler_notrimmer.supports #supports whatever the underlying Velvet binding supports

def assemble(path, accession, ref_url=None, k_value="177", fasta_file_name=None, callback=lambda s: None, update_callback=lambda d: None, round_1_range=100, round_1_step=20, round_2_range=14, round_2_step=2, debug=True, **parameters):
	
	callback("Velvet Optimize beginning: {} to {}".format(max(int(k_value) - round_1_range, 15), min(int(k_value) + round_1_range + round_1_step, 251)))
	start_time = datetime.datetime.today()
	
	tmp_path = tempfile.mkdtemp()
	if not fasta_file_name:
		fasta_file_name = "{}.optimized.fasta".format(accession)
		
	parameters['debug'] = debug
	
	try:
		results = list()
		for k in range(max(int(k_value) - round_1_range, 15), min(int(k_value) + round_1_range + round_1_step, 251), round_1_step):
			parameters['path'] = os.path.join(tmp_path, str(k))
			if not os.path.exists(parameters['path']):
				os.mkdir(parameters['path'])
			def callback_wrapper(s): #closure to wrap callback method, if any, from assembly_dispatch
				callback("Velvet Optimize round 1, k={} ({} to {}): {}".format(k, max(int(k_value) - round_1_range, 15), min(int(k_value) + round_1_range + round_1_step, 251), s))
			try:
				results.append((k, velvet_assembler_notrimmer.assemble(accession, k_value = k, fasta_file_name = "{}.k{}.fasta".format(accession, k), callback=callback_wrapper, **parameters), ))
			except ValueError as e:
				print e
# 		try:
# 			results.append(("WORST", worst_assembler.assemble(path=os.path.join(tmp_path, "WORST"), callback=callback, **parameters)))
# 		except Exception:
# 			callback("WORST assembler did not complete.")
		
		try:		
			[r[1].__setitem__('contiguity', (int(r[1]['n50']) * int(r[1]['num_bases']) / int(r[1]['num_contigs']))) for r in results]
		except IndexError as e:
			print "Index error on contiguity calculation - any results?"
			print results
		try:
			best = results[0][1]
			best_k = max(int(k_value) - round_1_range, 15)
		except IndexError:
			raise ValueError("Empty results list; no k_value produced an assembly.")
			
			

		for (k, result) in results:
			if result['contiguity'] > best['contiguity']:
				best = result
				best_k = k



		callback("Best k so far was {}; starting round 2".format(best_k))

		for k in range(max(best_k - round_2_range, 15), min(best_k + round_2_range + round_2_step, 251), round_2_step):
			parameters['path'] = os.path.join(tmp_path, str(k))
			if not os.path.exists(parameters['path']):
				os.mkdir(parameters['path'])
			def callback_wrapper(s): #closure to wrap callback method, if any, from assembly_dispatch
				callback("Velvet Optimize round 2, k={} ({} to {}): {}".format(k, max(int(best_k) - round_2_range, 15), min(int(best_k) + round_2_range + round_2_step, 251), s))
			try:
				results.append((k, velvet_assembler_notrimmer.assemble(accession, k_value = k, fasta_file_name = "{}.k{}.fasta".format(accession, k), callback=callback_wrapper, **parameters), ))
			except ValueError as e:
				print e
				
		try:		
			[r[1].__setitem__('contiguity', (int(r[1]['n50']) * int(r[1]['num_bases']) / int(r[1]['num_contigs']))) for r in results]
		except IndexError as e:
			print "Index error on contiguity calculation - any results?"
			print results
		#print results
		
		old_best_k = best_k
		
		for (k, result) in results:
			result['_k'] = k
			if result['contiguity'] > best['contiguity']:
				best = result
				best_k = k
				
		if best_k == max(old_best_k - round_2_range, 15) and not best_k == 15:
			callback("k-range wasn't broad enough; re-center and recursively try again.")
			def ucb(d):
				best_k = d['k_value']
			try:
				parameters['path'] = path
				best = assemble(accession, best_k, fasta_file_name, callback, ucb, round_1_range, round_1_step, round_2_range, round_2_step, **parameters)
			except IOError:
				print "Recursive subassembly terminated; out of room. {} will have to do.".format(best_k)
			except Exception:
				if debug:
					import traceback
					class CallbackWrapper(object):
						"duck-typed object to route stack traceback through callback system"
						def write(s):
							for line in s.split("\n"):
								callback(line)
					traceback.print_exc(CallbackWrapper())
		
		update_callback({'k_value':best_k, 'job_type':'Velvet'})
		shutil.copyfile(os.path.join(tmp_path, str(best_k), best['fasta_file']), os.path.join(path, fasta_file_name))
		
		import csv
		with open(os.path.join(path, "velvet_optimize.stats.txt"), 'w') as stats_file:
				wr = csv.DictWriter(stats_file, fieldnames=sorted(results[0][1].keys()), delimiter="\t")
				try:
					wr.writeheader()
				except TypeError:
					stats_file.write("\t".join(results[0][1].keys()))
				for (k, result) in sorted(results, key=lambda r: r[0]):
					try:
						wr.writerow(result)
					except TypeError:
						stats_file.write("\t".join(result.values()))
		best['elapsed'] = str(datetime.datetime.today() - start_time)
		callback("Velvet Optimize Winner: Kmer size {_k} with n50 {n50} and {num_contigs} contigs in {elapsed}".format(**best))
		try:
			results.sort(key=lambda r: r[0])
			for r in results:
				r[1]['fasta_file'] = "{}/{}".format(r[0], r[1]['fasta_file'])
				r[1]['assembler'] = "Velvet_{:03}".format(r[0])	
							
			quast_results = fasta_statter.quast_compare(tmp_path, [r[1] for r in results], callback=callback, update_callback=update_callback, debug=debug)
			for r in results:
				r[1].update(filter(lambda qr: r[1]['assembler'] in qr['assembler'], quast_results)[0])
			results = [r[1] for r in results]
			results.sort(key=lambda r: r['assembler'])
			try:
				callback("Copying {}...".format(os.path.join(tmp_path, 'quast/')))
				if os.path.exists(os.path.join(path, 'quast/')):
					shutil.rmtree(os.path.join(path, 'quast/'))
				shutil.copytree(os.path.join(tmp_path, 'quast/'), os.path.join(path, 'quast/'))
			except shutil.Error:
				pass
			except OSError:
				callback("QUAST copy failed.")
				if debug:
					import traceback
					import sys
					traceback.print_exc(sys.stdout)
			
			if debug:
				#headers = sorted(results[0].keys())
				headers = ('assembler',
						   'n50',
						   'num_contigs',
						   'Total length',
						   'Reference length',
						   'Genome fraction (%)',
						   '# misassemblies',
							'goodness')
				widths = [max([max(len(str(r.get(h, ''))), len(str(h))) for r in results]) + 2 for h in headers] #get longest string in each field in results
				callback("".join([str(h).ljust(w) for (w, h) in zip(widths, headers)]))
				for r in results:		
					callback("".join([str(r.get(h, '-')).ljust(w) for (w, h) in zip(widths, headers)]))
		except Exception:
			print "QUAST failed."
			if debug:
				import traceback
				class CallbackWrapper(object):
					"duck-typed object to route stack traceback through callback system"
					def write(s):
						for line in s.split("\n"):
							callback(line)
				traceback.print_exc(CallbackWrapper())
	finally:
		callback("Cleaning scratch space: {}".format(tmp_path))
		shutil.rmtree(tmp_path)

	
	best['fasta_file'] = fasta_file_name
	return best
	
	
if __name__ == "__main__":
	#debug
	import datetime
	def cb(d):
		print "[{}] {}".format(datetime.datetime.today().ctime(), d)
	def bcb(d):
		for (k, v) in d.items():
			cb("{} : {}".format(k, v))
	bcb(assemble(path='/shared/gn2/CFSANgenomes/CFSAN001659/asm/', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R1_001.fastq', 
			 reads2='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R2_001.fastq', 
			 accession='CFSAN001659_01', 
			 callback=cb,
			 update_callback=bcb,
			 k_value=144,
			 debug=True))
	bcb(assemble(path='/shared/gn2/CFSANgenomes/CFSAN001659/asm/', 
			 reads1='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R1_001.fastq', 
			 #reads2='/shared/gn2/CFSANgenomes/CFSAN001659/CFSAN001659_01/CFSAN001659_S11_L001_R2_001.fastq', 
			 accession='CFSAN001659_01', 
			 data_type='IonTorrent',
			 callback=cb,
			 update_callback=bcb,
			 k_value=144,
			 debug=True))