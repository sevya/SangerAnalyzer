#! /usr/bin/env python3

import PySimpleGUI as sg
from SangerAnalyzer import SangerAnalyzer

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

	## Remove window from the screen
	window.close()

	## Now do my heavy lifting
	# try:
	analyzer = SangerAnalyzer( out_path, parent_sequence, in_files=in_file_list, max_N=max_N, prob_limit=limit, width=60 )
	# sg.popup( 'Success!' )

	# except Exception as e:
	# 	sg.popup( 'Error: %s'%(e) )
	# 	exit()
