#! /usr/bin/env python3

from Bio import SeqIO, Seq
import numpy as np
import pandas as pd
import subprocess, os

def _phred_to_prob( Q ):
	return 10**(float(Q)/-10)

def _abi_to_seqs( in_files ):
	in_sequences = []
	for file in in_files:
		read = SeqIO.parse( file, 'abi' )
		name = os.path.split( file )[ 1 ].split('.')[0]

		for seq in read:
			in_sequences.append(( \
				name,             \
				str(seq.seq),     \
				seq.letter_annotations['phred_quality']) \
				)

	return in_sequences 

## CLC Workbench trimming algorithm as described here:
## http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/900/index.php?manual=Quality_trimming.html
def _trim_sequence( nt_seq, qual_arr, limit=0.05, max_N=2 ):
	run_sum = [0]

	limit_minus_p = [ limit - _phred_to_prob( base ) for base in qual_arr ]

	for value in limit_minus_p:
		run_sum.append( run_sum[-1] + value )
		if run_sum[-1] < 0:
			run_sum[-1] = 0

	run_sum = run_sum[1:]

	## if run_sum is all 0 (no good bases), return empty
	if run_sum == [0]*len(run_sum):
		return ('',[])

	start_idx = 0
	end_idx = np.argmax( run_sum )

	for ii in range( end_idx, 0, -1 ):
		if run_sum[ii] == 0:
			start_idx = ii+1
			break

	seq_trim = nt_seq[ start_idx:end_idx+1 ]
	qual_trim = qual_arr[ start_idx:end_idx+1 ]
	n_idx = [idx for idx in range(len(seq_trim)) if seq_trim[idx]=='N']

	if len(n_idx) > max_N:

		candidate_substrings = []
		idx_list_to_scan = [-1]+n_idx+[len(seq_trim)]
		for ii, idx1 in enumerate(idx_list_to_scan):
			for jj, idx2 in enumerate(idx_list_to_scan[ii+1:]):
				substr = seq_trim[idx1+1:idx2]
				N_count = sum([1 for ch in substr if ch=='N'])
				if N_count <= max_N:
					candidate_substrings.append( (substr,N_count,len(substr),idx1,idx2))
		
		candidate_substrings.sort( key=lambda a: a[2], reverse=True )

		trim_seq,N_count,length,idx1,idx2 = candidate_substrings[0]

		qual_trim = qual_trim[idx1+1:idx2]
		seq_trim = trim_seq

	return (''.join(seq_trim),qual_trim)

def _format_html( wt_seq, var_seq, primer, start_idx ):

	wt_array  = []
	var_array = []
	all_muts = []

	pos_it = start_idx
	for wt_nt, var_nt in zip(wt_seq,var_seq):
		wt_array.append( wt_nt )
		if wt_nt != var_nt:
			char = '<span style="background-color:red;font-weight:bold">%s</span>'%(var_nt)
			all_muts.append( '%s%d%s'%( wt_nt, pos_it, var_nt))
		else:
			char = var_nt

		var_array.append( char )
		## if WT NT is '-', that means there's an insertion, and I don't increment pos_it
		if wt_nt != '-':
			pos_it += 1


	wt_str  = 'WT'.ljust(20) + str(start_idx).rjust(4) + '  ' + ''.join( wt_array ) + '  ' + str(pos_it-1)
	var_str = primer.ljust(20) + str(start_idx).rjust(4) + '  ' + ''.join( var_array ) + '  ' + str(pos_it-1)

	html_str = '<pre>' + wt_str + '</br>' + var_str + '</br></br>' + '</pre>'
	return (html_str,all_muts)


## start_idx is zero-indexed
def _translate( seq, start_idx ):
	if start_idx%3 == 0:
		return str(Seq.Seq( seq.replace('-','').replace('_','') ).translate())
	elif start_idx%3 == 1:
		return str(Seq.Seq( seq.replace('-','').replace('_','')[2:] ).translate())
	elif start_idx%3 == 2:
		return str(Seq.Seq( seq.replace('-','').replace('_','')[1:] ).translate())
	else:
		return ''


# Yield successive chunks of size width from array.
def _divide_chunks( array, width ):
	# looping till length l
	for i in range(0, len(array), width): 
		yield array[i:i + width]

class SangerAnalyzer:


	## Constructor takes the following arguments:
	## in_files: a list of abi files to process
	## in_sequences: a list of tuples in the format (name,NT_sequence,quality_string) to process
	## One and only one of in_files and in_sequences should be provided
	## out_path: where to write output files
	## max_N: how many N nucleotides can be included in trimmed sequence
	## prob_limit: probability limit for trimming
	## width: width of lines in output html
	def __init__( self, out_path, parent_sequence, in_files=[], in_sequences=[], max_N=0, prob_limit=0.001, width=60 ):
		
		## Verify that one and only one of in_files and in_sequences is provided
		if not in_files and not in_sequences:
			raise ValueError('One and only one of in_files and in_sequences should be provided')
		elif len( in_files ) and len( in_sequences ):
			raise ValueError('One and only one of in_files and in_sequences should be provided')

		## Files are provided - turn them into sequences
		elif in_files != []:
			in_sequences = _abi_to_seqs( in_files )

		## ============================================================================ 
		## Ready to process data
		## My sequences and quality data are now in in_sequences list
		## I can go through and do the processing
		## ============================================================================ 



		## ============================================================================ 
		## Step 1. Make a blast db for the parent sequence, both NT and AA 
		## ============================================================================ 

		
		## Find executable based on location of script
		wd = os.path.dirname( __file__ )
		makeblastdb_exe = os.path.join( wd, 'blast-2.11.0+', 'bin', 'makeblastdb.exe' )
		## Make the commands 
		makeblastdb_nt_cmd = [ makeblastdb_exe, '-dbtype', 'nucl', '-parse_seqids', '-in', 'parent_nt.fasta', '-input_type', 'fasta', '-out', 'parent_nt' ]
		makeblastdb_aa_cmd = [ makeblastdb_exe, '-dbtype', 'prot', '-parse_seqids', '-in', 'parent_aa.fasta', '-input_type', 'fasta', '-out', 'parent_aa' ]
		## Write the parent sequences to fasta
		with open( 'parent_nt.fasta', 'w' ) as out:
			out.write( '>parent\n%s\n'%(parent_sequence) )
		parent_aa_sequence = _translate( parent_sequence, 0 )
		with open( 'parent_aa.fasta', 'w' ) as out:
			out.write( '>parent\n%s\n'%(parent_aa_sequence) )

		## Make the blast dbs
		subprocess.call( makeblastdb_nt_cmd )
		subprocess.call( makeblastdb_aa_cmd )

		## ============================================================================ 
		## Step 2. Trim the input NT sequences and write them to a fasta 
		## ============================================================================ 
		with open( 'sanger_nt.fasta', 'w' ) as out:
			for name, seq, qual in in_sequences:
				seq_trim, qual_trim = _trim_sequence( \
					seq, qual, \
					limit=prob_limit, \
					max_N=max_N )

				out.write( '>%s\n%s\n'%(name, seq_trim) )

		## ============================================================================ 
		## Step 3.Run BLAST of NT sequences against parent
		## ============================================================================ 
		blast_exe = os.path.join( wd, 'blast-2.11.0+', 'bin', 'blastn.exe' )
		
		blast_cmd = [ blast_exe, '-db', 'parent_nt', \
				'-out', 'sanger_nt.fasta.out', \
				'-query', 'sanger_nt.fasta', \
				'-outfmt', '6 qseqid qframe sframe length pident nident qstart qend sstart send qseq sseq sstrand'
				]

		## Remove previous out file
		if os.path.exists( 'sanger_nt.fasta.out' ):
		 	os.remove( 'sanger_nt.fasta.out' )
		
		## Set BLASTDB env variable
		os.environ['BLASTDB'] = wd
		## Run blast
		subprocess.call( blast_cmd )

		## ============================================================================ 
		## Step 4. Process the aligned sequences using a Pandas dataframe
		## ============================================================================ 
		df = pd.read_table( 'sanger_nt.fasta.out', header=None, \
			names=['name','qframe','sframe','length','pident','nident','qstart','qend','sstart','send','qseq','sseq','sstrand'])

		## TODO fix hack: i need to separate samples from primers by dashes. 
		## Replace the dash in T7-term. hope no other primers have this problem
		df['name'] = df['name'].apply( lambda a: a.replace('T7-Term','T7Term'))
		df['prefix'] = df['name'].apply( lambda a: '-'.join(a.split('-')[:-1]))

		## Figure out what are the different samples grouped by different primers
		samples = list(df['name'])
		## Separate the primer from the end of the name, then make each one unique
		prefixes = list(set([ '-'.join(s.split('-')[:-1]) for s in samples ]))
		primers = list(set([ s.split('-')[-1] for s in samples ]))
		prefixes.sort()
		primers.sort()

		## Go through each sample, grouped by prefixes, and see what the different primers look like
		aa_fasta_str = '' ## get a string ready to add in all the AA sequences in fasta format
		for prefix in prefixes:
			subdf = df[df['prefix']==( prefix )]
			all_muts = []
			all_html = []
			html_header = '<pre><header style="font-weight:bold">%s</header>'%(prefix)
			
			## Don't iterate by primers (all samples may not have each primer), just go through rows in subdf
			for ii, row in subdf.iterrows():
				primer_html = '' ## primer html holds the html string for each primer. 
							     ## All html is a list of all of the different primer html strings
				primer = row['name'].split('-')[-1]
				## If this primer is a reverse primer, reverse complement the sequence
				## Now everything will be in forward orientation
				if row['sstrand'] == 'minus':
					qseq = str(Seq.Seq( row['qseq'] ).reverse_complement())
					sseq = str(Seq.Seq( row['sseq'] ).reverse_complement())
					sstart = int(row['send'])
				else:
					qseq = row['qseq']
					sseq = row['sseq']
					sstart = row['sstart']

				## Positional iterator is where am I in WT sequence - given by sstart
				start_idx = sstart
				## break sequence into chunks of preset width to format as html
				## qseq is variant, sseq is wt
				wt_chunks = list(_divide_chunks( sseq, width ))
				var_chunks = list(_divide_chunks( qseq, width ))

				for chunk1,chunk2 in zip(wt_chunks,var_chunks):
					html,muts = _format_html( chunk1, chunk2, primer, start_idx )
					start_idx += width
					all_muts += muts
					primer_html += html

				all_html.append( (sstart, primer_html ))

				## write AA sequences to fasta file
				aa_seq = _translate( qseq, sstart-1 )
				aa_fasta_str += '>%s\n%s\n'%(row['name'],aa_seq)


			## sort and deduplicate all muts
			all_muts = list(set(all_muts))
			all_muts.sort( key=lambda a: int(a[1:-1]))

			## Sort the html for each primer by where the start index is
			all_html.sort( key=lambda a: a[0] )
			html_footer = '<header style="font-weight:bold">All mutations</header><span style="color:red;font-weight:bold">%s</span></pre>'%(';'.join(all_muts))

			## demarcate the region between primers
			demarc = '='*(width + 30) + '</br></br>'

			html_combined = html_header + demarc.join( [tup[1] for tup in all_html] ) + demarc + html_footer

			out_file = os.path.join( out_path, prefix+'_NT.html' )
			with open( out_file,'w') as out:
				out.write( html_combined )

		

			
		## ============================================================================ 
		## Step 5. Blast the protein sequences
		## ============================================================================ 
		## Write out fasta file with AA sequences
		with open( 'sanger_aa.fasta','w') as out:
			out.write( aa_fasta_str )

		blast_exe = os.path.join( wd, 'blast-2.11.0+', 'bin', 'blastp.exe' )
		
		blast_cmd = [ blast_exe, '-db', 'parent_aa', \
				'-out', 'sanger_aa.fasta.out', \
				'-query', 'sanger_aa.fasta', \
				'-max_hsps', '1', \
				'-outfmt', '6 qseqid qframe sframe length pident nident qstart qend sstart send qseq sseq'
				]

		## Remove previous out file
		if os.path.exists( 'sanger_aa.fasta.out' ):
		 	os.remove( 'sanger_aa.fasta.out' )
		
		## Set BLASTDB env variable
		os.environ['BLASTDB'] = wd
		## Run blast
		subprocess.call( blast_cmd )


		## ============================================================================ 
		## Step 6. Process the aligned aa sequences using a Pandas dataframe
		## ============================================================================ 
		df = pd.read_table( 'sanger_aa.fasta.out', header=None, \
			names=['name','qframe','sframe','length','pident','nident','qstart','qend','sstart','send','qseq','sseq'])

		## hack: i need to separate samples from primers by dashes. Replace the dash in T7-term. hope no other primers have this problem
		df['name'] = df['name'].apply( lambda a: a.replace('T7-Term','T7Term'))
		df['prefix'] = df['name'].apply( lambda a: '-'.join(a.split('-')[:-1]))

		for prefix in prefixes:
			subdf = df[df['prefix']==prefix]

			all_muts = []
			all_html = [] 
			html_header = '<pre><header style="font-weight:bold">%s</header>'%(prefix)
			
			

			## Don't iterate by primers (all samples may not have each primer), just go through rows in subdf
			for ii, row in subdf.iterrows():
				
				primer = row['name'].split('-')[-1]
				primer_html = '' ## primer html holds the html string for each primer. 
							     ## All html is a list of all of the different primer html strings

				## Positional iterator is where am I in WT sequence - given by sstart
				start_idx = row['sstart']
				## break sequence into chunks of preset width to format as html
				## qseq is variant, sseq is wt
				wt_chunks = list(_divide_chunks( row['sseq'], width ))
				var_chunks = list(_divide_chunks( row['qseq'], width ))

				for chunk1,chunk2 in zip(wt_chunks,var_chunks):
					html,muts = _format_html( chunk1, chunk2, primer, start_idx )
					start_idx += width
					all_muts += muts
					primer_html += html


				all_html.append( (row['sstart'], primer_html ))
				primer_html = ''


			## sort and deduplicate all muts
			all_muts = list(set(all_muts))
			all_muts.sort( key=lambda a: int(a[1:-1]))

			## Sort the html for each primer by where the start index is
			all_html.sort( key=lambda a: a[0] )
			html_footer = '<header style="font-weight:bold">All mutations</header><span style="color:red;font-weight:bold">%s</span></pre>'%(';'.join(all_muts))

			## demarcate the region between primers
			demarc = '='*(width + 30) + '</br></br>'

			html_combined = html_header + demarc.join( [tup[1] for tup in all_html] ) + demarc + html_footer

			## Write it out
			out_file = os.path.join( out_path, prefix+'_AA.html' )
			with open( out_file,'w') as out:
				out.write( html_combined )


