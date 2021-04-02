#! /usr/bin/env python3

from Bio import SeqIO, Seq
import numpy as np
import pandas as pd
import subprocess, os
import PySimpleGUI as sg

def smooth( array, window=5 ):
	ret_array = [0]*int((window-1)/2)
	
	for ii in range( 0, len(array)-window ):
		ret_array.append( np.mean( array[ii:ii+window] ))

	ret_array += [0]*int((window-1)/2)
	return ret_array

def phred_to_prob( Q ):
	return 10**(Q/-10)

## CLC Workbench trimming algorithm as described here:
## http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/900/index.php?manual=Quality_trimming.html
def trim_sequence( nt_seq, qual_arr, limit=0.05, max_N=2 ):
	run_sum = [0]

	limit_minus_p = [ limit - phred_to_prob( base ) for base in qual_arr ]

	for value in limit_minus_p:
		run_sum.append( run_sum[-1] + value )
		if run_sum[-1] < 0:
			run_sum[-1] = 0

	run_sum = run_sum[1:]

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

def format_html( wt_seq, var_seq, primer_name, start_idx ):

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



def translate( seq, start_idx ):
	if start_idx%3 == 0:
		return str(Seq.Seq( seq.replace('-','').replace('_','') ).translate())
	elif start_idx%3 == 1:
		return str(Seq.Seq( seq.replace('-','').replace('_','')[2:] ).translate())
	elif start_idx%3 == 2:
		return str(Seq.Seq( seq.replace('-','').replace('_','')[1:] ).translate())
	else:
		return ''


# Yield successive chunks of size width from array.
def divide_chunks( array, width ):
	# looping till length l
	for i in range(0, len(array), width): 
		yield array[i:i + width]
  

def is_float( string ):
	try:
		float(string)
		return True
	except:
		return False


def is_int( string ):
	try:
		int(string)
		return True
	except:
		return False

if __name__=='__main__':

	layout = [  [sg.Text('Input NT sequence of parent:')],
				[sg.Multiline('',size=(45,5))],
				[sg.Text('Select sequence files in .abi format'), sg.FilesBrowse(file_types=(('AB1 files','*.ab1'),))],
				[sg.Text('Choose where to save output files'), sg.FolderBrowse()],
				[sg.Text('Trimming options')],
				[sg.Text('Maximum N allowed in sequence'), sg.InputText('0')],
				[sg.Text('Minimum quality allowed'), sg.InputText('0.001')],
				[sg.Submit(), sg.Cancel() ]
			]

	## Create the window
	window = sg.Window( 'Sanger sequence analyzer', layout )

	## Create a while loop to keep showing window until all params are satisfied
	all_params_input = False
	
	while not all_params_input:
		# Display and interact with the Window
		event, values = window.read()

		parent_sequence = values[ 0 ]
		in_files = values[ 'Browse' ]
		out_path = values[ 'Browse0' ]
		max_N = values[ 1 ]
		limit = values[ 2 ]


		all_params_input = True
		## Check that all parameters were input
		if parent_sequence.strip() == '':
			sg.popup( 'Please enter a parent sequence' )
			all_params_input = False
			continue

		if in_files.strip() == '':
			sg.popup( 'Please select AB1 files to process' )
			all_params_input = False
			continue

		if out_path.strip() == '':
			sg.popup( 'Please choose where to save the output files' )
			all_params_input = False
			continue

		## Check that max_N and limit are ints and floats (if left alone they will be default values)
		if not is_int( max_N ):
			sg.popup( 'Maximum N should be an integer' )
			all_params_input = False
			continue

		if not is_float( limit ):
			sg.popup( 'Minimum quality should be a numeric value' )
			all_params_input = False
			continue



	in_file_list = in_files.split(';') ## TODO uncomment this line - just for testing purposes
	max_N = int( max_N )
	limit = float( limit )
	

	limit = 0.001
	max_N = 0
	# ## For testing purposes - don't invoke GUI to get these values. TODO comment them out
	# in_file_list = 'C:/Users/Alex/OneDrive - Singular Genomics Systems/Desktop/Souad sequencing data/c8-T7-Term.ab1'.split(';')

	# out_path = 'C:/Users/Alex/OneDrive - Singular Genomics Systems/Desktop/Souad sequencing data/'
	# parent_sequence = 'ATGATTCTGGACGCTGATTATATTACTGAAGATGGTAAACCGATTATTCGTATTTTTAAAAAAGAAAATGGCGAGTTCAAAGTTGAATATGACCGTAACTTTCGTCCGTACATCTACGCGCTGTTGCGCGACGATAGCGCGATCGATGAGATTAAGAAAATTACCGCGCAGCGTCATGGTAAAGTTGTTCGCATCGTTGAAACCGAGAAAATTCAACGTAAATTCCTGGGCCGCCCAATTGAAGTGTGGAAGCTGTACCTGGAGCATCCGCAAGATGTCCCGGCGATCCGTGACAAGATTCGCGAGCACCCGGCCGTCGTCGACATTTTCGAATACGATATTCCGTTCGCAAAGCGTTACCTGATCGATAAGGGTCTGACCCCGgcaGAGGGTAATGAAAAGCTGACGTTCCTGGCTGTCgcaATTgcggcgTTGTACCACGAGGGTGAAGAGTTTGGTAAGGGCCCGGTCATTATGATCAGCTACGCGGATGAAGAGGGCGCCAAAGTTATCACGTGGAAAAAAATTGATCTGCCGTACGTTGAAGTTGTGTCCAGCGAGCGCGAGATGATTAAACGCTTGATTCGTGTGATTAAAGAAAAAGATCCAGACGTGATCATTACCTATAATGGTGACAACTTTGACTTTCCGTACTTGCTGAAACGTGCTGAGAAACTGGGTATCAAGCTGTTGCTGGGTCGCGATAATAGCGAGCCGAAGATGCAAAAAATGGGCGATAGCCTGGCAGTCGAGATCAAGGGTCGCATCCACTTTGATCTCTTTCCGGTGATTCGTCGCACGATCAATCTGCCGACCTATACGCTGGAAGCTGTCTACGAGGCAATCTTTGGTAAGCCGAAAGAAAAAGTCTATGCGGACGAAATTGCGAAAGCGTGGGAAACCGGCGAGGGCCTGGAGCGTGTGGCAAAGTACTCTATGGAAGATGCCAAAGTGACCTATGAACTGGGTCGTGAGTTCTTCCCAATGGAAGCCCAGTTGGCGCGCTTGGTGGGCCAACCGGTTTGGGACGTTTCCCGTAGCAGCACCGGTAACCTGGTTGAGTGGTTTCTGTTGCGTAAAGCGTATGAGCGTAATGAACTGGCACCGAACAAGCCTGACGAGAAAGAATATGAACGTCGCCTGCGTGAATCTTACGAGGGTGGTTACGTCAAAGAACCGGAAAAGGGTCTGTGGGAAGGCATCGTGAGCCTGGATTTCCGTAGCgcaggcccgAGCATCATCATCACGCACAATGTTAGCCCGGACACCCTGAACCGCGAGGGCTGCGAAGAGTACGACGTTGCGCCGAAAGTCGGCCATCGTTTTTGTAAAGACTTCCCTGGTTTCATCCCAAGCCTGCTGGGTCAGCTGCTGGAAGAGAGACAGAAAATTAAAAAACGCATGAAAGAATCGAAAGATCCGGTTGAGAAAAAGCTGCTGGATTACCGCCAGCGTgtgATCAAGATTCTGGCTAACTCATATTATGGCTACTACGGTTATGCTAAAGCGCGTTGGTACTGTAAAGAGTGCGCGGAGTCCGTCtctGCGTGGGGTCGCCAGTATATCGATCTGGTGCGTCGCGAGCTGGAAGCGCGTGGTTTTAAGGTCCTGTACATCGATACTGACGGTCTGTATGCAACCATCCCTGGTGTCAAAGACTGGGAAGAGGTTAAGCGTCGTGCACTGGAATTTGTGGACTATATCAATTCTAAGTTGCCGGGTGTGCTGGAGCTGGAGTACGAAGGCTTCTATGCACGCGGCTTTTTCGTTacgAAAAAGAAATACGCACTGATCGACGAAGAGGGCAAGATTGTGACTCGTGGTCTGGAAATCGTTCGTCGCGACTGGAGCGAGATTGCAAAAGAAACCCAAGCTCGCGTTCTGGAAGCAATCCTGAAACATGGTAACGTCGAAGAActgGTCAAGATCGTGAAAGATGTCACCGAAAAGTTGACCAACTACGAAGTTCCACCGGAAAAACTGGTGATTTATGAGCAAATCACGCGTCCGATCAATGAATATAAGGCCATTGGCCCGCACGTCGCGGTGGCCAAGCGCCTGATGGCGCGTGGTATCAAAGTGAAACCGGGTATGGTTATTGGTTACATCGTGCTGCGTGGCGACGGCCCGATTAGCAAACGTGCGATCAGCATTGAAGAATTTGACCCGCGTAAGCACAAATATGACGCGGAATACTATATCGAGAATCAAGTGCTGCCGGCCGTGGAACGCATTCTGAAAGCTTTCGGCTACAAGCGTGAAGATTTGCGCTGGCAGAAAACCAAACAGGTTGGTCTTGGTGCGTGGATCAAGGTCAAAAAGTCCTAA'
	## Remove window from the screen
	window.close()

	## Now do my heavy lifting

	try:
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
		parent_aa_sequence = translate( parent_sequence, 0 )
		with open( 'parent_aa.fasta', 'w' ) as out:
			out.write( '>parent\n%s\n'%(parent_aa_sequence) )

		## Make the blast dbs
		subprocess.call( makeblastdb_nt_cmd )
		subprocess.call( makeblastdb_aa_cmd )

		## ============================================================================ 
		## Step 2. Trim the input NT sequences and write them to a fasta 
		## ============================================================================ 

		with open( 'sanger_nt.fasta', 'w' ) as out:
			for file in in_file_list:
				read = SeqIO.parse( file, 'abi' )
				name = os.path.split( file )[ 1 ].split('.')[0]

				for seq in read:

					seq_trim, qual_trim = trim_sequence( \
						str(seq.seq), \
						seq.letter_annotations['phred_quality'], \
						limit=limit, \
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
		## hack: i need to separate samples from primers by dashes. Replace the dash in T7-term. hope no other primers have this problem
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

			## Instead of making contig map, just annotate my mutations as I go and add them to a list
			
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
				width = 60
				start_idx = sstart
				## break sequence into chunks of preset width to format as html
				## qseq is variant, sseq is wt
				wt_chunks = list(divide_chunks( sseq, width ))
				var_chunks = list(divide_chunks( qseq, width ))

				for chunk1,chunk2 in zip(wt_chunks,var_chunks):
					html,muts = format_html( chunk1, chunk2, primer, start_idx )
					start_idx += width
					all_muts += muts
					primer_html += html

				all_html.append( (sstart, primer_html ))

				## write AA sequences to fasta file
				aa_seq = translate( qseq, sstart-1 )
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
				width = 60
				start_idx = row['sstart']
				## break sequence into chunks of preset width to format as html
				## qseq is variant, sseq is wt
				wt_chunks = list(divide_chunks( row['sseq'], width ))
				var_chunks = list(divide_chunks( row['qseq'], width ))

				for chunk1,chunk2 in zip(wt_chunks,var_chunks):
					html,muts = format_html( chunk1, chunk2, primer, start_idx )
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

		## ============================================================================ 
		## Done!!
		## ============================================================================ 
		sg.popup( 'Success!' )

	except Exception as e:
		sg.popup( 'Error: %s'%(e) )
		exit()
