#! /usr/bin/env python3

import SangerAnalyzer as sa
import unittest, os

class TestSangerAnalyzer( unittest.TestCase ):

	def test_phred_to_prob( self ):
		self.assertEqual( sa._phred_to_prob( 30 ), 0.001 )
		self.assertEqual( sa._phred_to_prob( 10 ), 0.1 )
		self.assertEqual( sa._phred_to_prob( 0 ), 1 )

		phreds = [4, 8, 15, 30, 32, 23, 10, 31, 31, 20, 28, 29, 30, 10, 31, 31, 20, 28, 32, 31, 32, 33, 29, 30, 31, 32, 32, 31, 32, 33, 33, 32, 33, 34, 30, 31, 32, 33, 33, 32, 33, 34, 23, 7, 7, 12]
		quals = [0.3981, 0.1585, 0.0316, 0.0010, 0.0006, 0.0050, 0.1000, 0.0008, 0.0008, 0.0100, 0.0016, 0.0013, 0.0010, 0.1000, 0.0008, 0.0008, 0.0100, 0.0016, 0.0006, 0.0008, 0.0006, 0.0005, 0.0013, 0.0010, 0.0008, 0.0006, 0.0006, 0.0008, 0.0006, 0.0005, 0.0005, 0.0006, 0.0005, 0.0004, 0.0010, 0.0008, 0.0006, 0.0005, 0.0005, 0.0006, 0.0005, 0.0004, 0.0050, 0.1995, 0.1995, 0.0631]

		for p,q in zip( phreds, quals ):
			self.assertAlmostEqual( sa._phred_to_prob( p ), q, places=4 )

	def test_trimming( self ):
		## ================== Regular trimming ==================
		sequence = 'GTGCACGTACTGGGTACTGGGTACTACGTACTGGGTATACTGGGTA'
		qual_arr = [4, 8, 15, 30, 32, 23, 10, 31, 31, 20, 28, 29, 30, 10, 31, 31, 20, 28, 32, 31, 32, 33, 29, 30, 31, 32, 32, 31, 32, 33, 33, 32, 33, 34, 30, 31, 32, 33, 33, 32, 33, 34, 23, 7, 7, 12]

		seq_trim, qual_trim = sa._trim_sequence( sequence, qual_arr, limit=0.05, max_N=2 )
		self.assertEqual( seq_trim, 'GCACGTACTGGGTACTGGGTACTACGTACTGGGTATACTGG' )
		self.assertEqual( qual_trim, \
			[15,30,32,23,10,31,31,20,28,29,30,10,31,31,20,28,32,31,32,33,29,30,31,32,32,31,32,33,33,32,33,34,30,31,32,33,33,32,33,34,23]
			)

		## ================== Regular trimming w alternate limit ==================
		sequence = 'GTGCACGTACTGGGTACTGGGTACTACGTACTGGGTATACTGGGTA'
		qual_arr = [4, 8, 15, 30, 32, 23, 10, 31, 31, 20, 28, 29, 30, 10, 31, 31, 20, 28, 32, 31, 32, 33, 29, 30, 31, 32, 32, 31, 32, 33, 33, 32, 33, 34, 30, 31, 32, 33, 33, 32, 33, 34, 23, 7, 7, 12]

		seq_trim, qual_trim = sa._trim_sequence( sequence, qual_arr, limit=0.001, max_N=2 )
		self.assertEqual( seq_trim, 'GGGTACTACGTACTGGGTATACTG' )
		self.assertEqual( qual_trim, \
			[32,31,32,33,29,30,31,32,32,31,32,33,33,32,33,34,30,31,32,33,33,32,33,34,]
			)

		## ================== N limit trimming - left side ==================
		sequence = 'GTGCACGTACTNGGTACTGGGTACTACGTACTGGGTATACTGGGTA'

		seq_trim, qual_trim = sa._trim_sequence( sequence, qual_arr, limit=0.05, max_N=0 )
		self.assertEqual( seq_trim, 'GGTACTGGGTACTACGTACTGGGTATACTGG' )
		self.assertEqual( qual_trim, \
			[30,10,31,31,20,28,32,31,32,33,29,30,31,32,32,31,32,33,33,32,33,34,30,31,32,33,33,32,33,34,23]
			)

		## ================== N limit trimming - right side ==================
		sequence = 'GTGCACGTACTGGGTACTGGGTACTACGTACTGGNTATACTGGGTA'
		           #  GCACGTACTGGGTACTGGGTACTACGTACTGGNTATACTGG              Sequence after quality trimming
		           #  GCACGTACTGGGTACTGGGTACTACGTACTGG                       Sequence after N trimming


		seq_trim, qual_trim = sa._trim_sequence( sequence, qual_arr, limit=0.05, max_N=0 )
		self.assertEqual( seq_trim, 'GCACGTACTGGGTACTGGGTACTACGTACTGG' )
		self.assertEqual( qual_trim, \
			[15,30,32,23,10,31,31,20,28,29,30,10,31,31,20,28,32,31,32,33,29,30,31,32,32,31,32,33,33,32,33,34]
			)

	    ## ================== N limit trimming - identify correct way to get longest read ==================
		sequence = 'GTGCACGTACTGGGTANTGGGTACTANGTACGNGGTATACTGGGTA'
		           #  GCACGTACTGGGTANTGGGTACTANGTACGNGGTATACTGG              Sequence after quality trimming
		           #  GCACGTACTGGGTANTGGGTACTA                               Solution 1 to N trimming (24 chars)
		           #                 TGGGTACTANGTACG                         Solution 2 to N trimming (15 chars)
		           #                           GTACGNGGTATACTGG              Solution 3 to N trimming (16 chars)

		seq_trim, qual_trim = sa._trim_sequence( sequence, qual_arr, limit=0.05, max_N=1 )
		self.assertEqual( seq_trim, 'GCACGTACTGGGTANTGGGTACTA' )
		self.assertEqual( qual_trim, \
			[15,30,32,23,10,31,31,20,28,29,30,10,31,31,20,28,32,31,32,33,29,30,31,32]
			)


		## ================== N limit trimming - identify correct way to get longest read (in the middle) ==================
		sequence = 'GTGCANGTACTGGGTANTGGNTACTANGTACGNGGTANACTGGGTA'
		           #  GCANGTACTGGGTANTGGNTACTANGTACGNGGTANACTGG              Sequence after quality trimming
		           #  GCANGTACTGGGTANTGG                                     Solution 1 to N trimming (18 chars)
		           #      GTACTGGGTANTGGNTACTA                               Solution 2 to N trimming (20 chars)
		           #                 TGGNTACTANGTACG                         Solution 3 to N trimming (15 chars)
		           #                     TACTANGTACGNGGTA                    Solution 4 to N trimming (16 chars)
		           #                           GTACGNGGTANACTGG              Solution 5 to N trimming (16 chars)

		seq_trim, qual_trim = sa._trim_sequence( sequence, qual_arr, limit=0.05, max_N=2 )
		self.assertEqual( seq_trim, 'GTACTGGGTANTGGNTACTA' )
		self.assertEqual( qual_trim, \
			[10,31,31,20,28,29,30,10,31,31,20,28,32,31,32,33,29,30,31,32]
			)

		## ================== One base left ==================
		sequence = 'GTGCACGTACTGGGTACTGGGTACTACGTACTGGGTATACTGGGTA'
		qual_arr = [4]*10 + [15] + [13]*35

		seq_trim, qual_trim = sa._trim_sequence( sequence, qual_arr, limit=0.05, max_N=2 )
		self.assertEqual( seq_trim, 'T' )
		self.assertEqual( qual_trim, [15] )

		## ================== No bases left ==================
		sequence = 'GTGCACGTACTGGGTACTGGGTACTACGTACTGGGTATACTGGGTA'
		qual_arr = [4]*10 + [12] + [13]*35

		seq_trim, qual_trim = sa._trim_sequence( sequence, qual_arr, limit=0.05, max_N=2 )
		self.assertEqual( seq_trim, '' )
		self.assertEqual( qual_trim, [] )

	def test_translate( self ):
		sequence = 'ATGAAACATATGCCGAGAAAGATGTATAGTTGTGACTTTGAGACAACTACTAAAGTGGAAGACTGTAGGGTATGGGCGTATGGTTATATGAATATAGAAGATCACAGTGAGTACAAAATAGGTAATAGCCTGGATGAGTTTATGGCGTGGGTGTTGAAGGTACAAGCTGATC'

		## in frame
		self.assertEqual( sa._translate( sequence, 0 ),     'MKHMPRKMYSCDFETTTKVEDCRVWAYGYMNIEDHSEYKIGNSLDEFMAWVLKVQAD' )
		self.assertEqual( sa._translate( sequence[3:], 3 ),  'KHMPRKMYSCDFETTTKVEDCRVWAYGYMNIEDHSEYKIGNSLDEFMAWVLKVQAD' )
		self.assertEqual( sa._translate( sequence[6:], 6 ),   'HMPRKMYSCDFETTTKVEDCRVWAYGYMNIEDHSEYKIGNSLDEFMAWVLKVQAD' )

		sequence = 'TGAAACATATGCCGAGAAAGATGTATAGTTGTGACTTTGAGACAACTACTAAAGTGGAAGACTGTAGGGTATGGGCGTATGGTTATATGAATATAGAAGATCACAGTGAGTACAAAATAGGTAATAGCCTGGATGAGTTTATGGCGTGGGTGTTGAAGGTACAAGCTGATC'

		## frame +1
		self.assertEqual( sa._translate( sequence, 1 ),     'KHMPRKMYSCDFETTTKVEDCRVWAYGYMNIEDHSEYKIGNSLDEFMAWVLKVQAD' )
		self.assertEqual( sa._translate( sequence[3:], 4 ),  'HMPRKMYSCDFETTTKVEDCRVWAYGYMNIEDHSEYKIGNSLDEFMAWVLKVQAD' )
		self.assertEqual( sa._translate( sequence[6:], 7 ),   'MPRKMYSCDFETTTKVEDCRVWAYGYMNIEDHSEYKIGNSLDEFMAWVLKVQAD' )

		sequence = 'GAAACATATGCCGAGAAAGATGTATAGTTGTGACTTTGAGACAACTACTAAAGTGGAAGACTGTAGGGTATGGGCGTATGGTTATATGAATATAGAAGATCACAGTGAGTACAAAATAGGTAATAGCCTGGATGAGTTTATGGCGTGGGTGTTGAAGGTACAAGCTGATC'

		## frame +2
		self.assertEqual( sa._translate( sequence, 2 ), 'KHMPRKMYSCDFETTTKVEDCRVWAYGYMNIEDHSEYKIGNSLDEFMAWVLKVQAD' )
		self.assertEqual( sa._translate( sequence[3:], 5 ),  'HMPRKMYSCDFETTTKVEDCRVWAYGYMNIEDHSEYKIGNSLDEFMAWVLKVQAD' )
		self.assertEqual( sa._translate( sequence[6:], 8 ),   'MPRKMYSCDFETTTKVEDCRVWAYGYMNIEDHSEYKIGNSLDEFMAWVLKVQAD' )

	def test_chunks( self ):
		array = [1,2,3,4,5,6,7,8,9,10,11,12]
		self.assertEqual( list(sa._divide_chunks( array, 2 )),  [[1,2],[3,4],[5,6],[7,8],[9,10],[11,12]] )
		self.assertEqual( list(sa._divide_chunks( array, 3 )),  [[1,2,3],[4,5,6],[7,8,9],[10,11,12]] )
		self.assertEqual( list(sa._divide_chunks( array, 4 )),  [[1,2,3,4],[5,6,7,8],[9,10,11,12]] )
		self.assertEqual( list(sa._divide_chunks( array, 5 )),  [[1,2,3,4,5],[6,7,8,9,10],[11,12]] )
		self.assertEqual( list(sa._divide_chunks( array, 10 )), [[1,2,3,4,5,6,7,8,9,10],[11,12]] )
		self.assertEqual( list(sa._divide_chunks( array, 13 )), [[1,2,3,4,5,6,7,8,9,10,11,12]] )

	def test_abi_parser( self ):

		output = sa._abi_to_seqs( ['resources\\HD-1-F0_303_F.ab1','resources\\HD-1-T7.ab1','resources\\HD-1-T7-Term.ab1','resources\\HD-2-F0_303_F.ab1','resources\\HD-2-T7.ab1','resources\\HD-2-T7-Term.ab1'])

		for (a1,b1,c1),(a2,b2,c2) in zip( self.template, output ):
			self.assertEqual( a1, a2 )
			self.assertEqual( b1, b2 )
			self.assertEqual( c1, c2 )


	def test_protocol( self ):

		out_files = [ 'resources\\HD-%d_%s.html'%(ii,aa) for ii in range(1,3) for aa in ('NT','AA') ]

		for file in out_files:
			if os.path.isfile( file ):
				os.remove( file )


		parent = 'ATGAAACATATGCCGAGAAAGATGTATAGTTGTGACTTTGAGACAACTACTAAAGTGGAAGACTGTAGGGTATGGGCGTATGGTTATATGAATATAGAAGATCACAGTGAGTACAAAATAGGTAATAGCCTGGATGAGTTTATGGCGTGGGTGTTGAAGGTACAAGCTGATCTATATTTCCATAACCTCAAATTTGACGGAGCTTTTATCATTAACTGGTTGGAACGTAATGGTTTTAAGTGGTCGGCTGACGGATTGCCAAACACATATAATACGATCATATCTCGCATGGGACAATGGTACATGATTGATATATGTTTAGGCTACAAAGGGAAACGTAAGATACATACAGTGATATATGACAGCTTAAAGAAACTACCGTTTCCTGTTAAGAAGATAGCTAAAGACTTTAAACTAACTGTTCTTAAAGGTGATATTGATTACCACAAAGAAAGACCAGTCGGCTATAAGATAACACCCGAAGAATACGCCTATATTAAAAACGATATTCAGATTATTGCGGAAGCTCTGTTAATTCAGTTTAAGCAAGGTTTAGACCGGATGACAGCAGGCAGTGACAGTCTAAAAGGTTTCAAGGATATTATAACCACTAAGAAATTCAAAAAGGTGTTTCCTACATTGAGTCTTGGACTCGATAAGGAAGTGAGATACGCCTATAGAGGTGGTTTTACATGGTTAAATGATAGGTTCAAAGAAAAAGAAATCGGAGAAGGCATGGTCTTCGATGTTAATAGTCTATATCCTGCACAGATGTATAGCCGTCTCCTTCCATATGGTGAACCTATAGTATTCGAGGGTAAATACGTTTGGGACGAAGATTACCCACTACACATACAGCATATCAGATGTGAGTTCGAATTGAAAGAGGGCTATATACCCACTATACAGATAAAAAGAAGTAGGTTTTATAAAGGTAATGAGTACCTAAAAAGTAGCGGCGGGGAGATAGCCGACCTCTGGTTGTCAAATGTAGACCTAGAATTAATGAAAGAACACTACGATTTATATAACGTTGAATATATCAGCGGCTTAAAATTTAAAGCAACTACAGGTTTGTTTAAAGATTTTATAGATAAATGGACGTACATCAAGACGACATCAGAAGGAGCGATCAAGCAACTAGCAAAACTGATGTTAAACAGTCTATACGGTAAATTCGCTAGTAACCCTGATGTTACAGGGAAAGTCCCTTATTTAAAAGAGAATGGGGCGCTAGGTTTCAGACTTGGAGAAGAGGAAACAAAAGACCCTGTTTATACACCTATGGGCGTTTTCATCACTGCATGGGCTAGATACACGACAATTACAGCGGCACAGGCTTGTTATGATCGGATAATATACTGTGATACTGACAGCATACATTTAACGGGTACAGAGATACCTGATGTAATAAAAGATATAGTTGACCCTAAGAAATTGGGATACTGGGCACATGAAAGTACATTCAAAAGAGCTAAATATCTGAGACAGAAGACCTATATACAAGACATCTATATGAAAGAAGTAGATGGTAAGTTAGTAGAAGGTAGTCCAGATGATTACACTGATATAAAATTTAGTGTTAAATGTGCGGGAATGACTGACAAGATTAAGAAAGAGGTTACGTTTGAGAATTTCAAAGTCGGATTCAGTCGGAAAATGAAGCCTAAGCCTGTGCAAGTGCCGGGCGGGGTGGTTCTGGTTGATGACACATTCACAATCAAAGGCGGTAGCGGTCATCATCATCACCACCACTAA'

		obj = sa.SangerAnalyzer( 'resources', parent, in_sequences=self.template, max_N=0, prob_limit=0.001, width=60 )

		for file in out_files:
			self.assertTrue( os.path.isfile( file ) )

	def setUp(self):
		self.template = correct = [\
		('HD-1-F0_303_F', \
			'NNNNNNNNNNNTNNNNTCANGGGAACGTAAGATACATACAGTGATATATGACAGCTTAAAGAAACTACCGTTTCCTGTTAAGAAGATAGCTAAAGACTTTAAACTAACTGTTCTTAAAGGTGATATTGATTACCACAAAGAAAGACCAGTCGGCTATAAGATAACACCCGAAGAATACGCCTATATTAAAAACGATATTCAGATTATTGCGGAAGCTCTGTTAATTCAGTTTAAGCAAGGTTTAGACCGGATGACAGCAGGCAGTGACAGTCTAAAAGGTTTCAAGGATATTATAACCACTAAGAAATTCAAAAAGGTGTTTCCTACATTGAGTCTTGGACTCGATAAGGAAGTGAGATACGCCTATAGAGGTGGTTTTACATGTTAAATGATAGGTTCAAAGAAAAAGAAATCGGAGAAGGCATGGTCTTCGATGTTAATAGTCTATATCCTGCACAGATGTATAGCCGTCTCCTTCCATATGGTGAACCTATAGTATTCGAGGGTAAATACGTTTGGGACGAAGATTACCCACTACACATACAACATATCAGATGTGAGTTCGAATTGAAAGAGGGCTATATACCCACTATACAGATAAAAAGAAGTAGGTTTTATAAAGGTAATGTGTACCTAAAAAGTAGCGGCGGGGAGATAGCCGACCTCTGGTTGTCAAATGTAGACCTAGAATTAATGAAAGAACACTACGATTTATATAACGTTGAATATATCAGCGGCTTAAAATTTAAAGCAACTACAGGTTTGTTTAAAGATTTTATAGATAAATGGACGTACATCATGACGACATCAGTAGGAGCGATCAAGCAACTAGCAAAACTGATGTTAAACAGTCTA', \
			[1,3,1,4,7,5,5,6,6,6,6,11,7,4,7,8,15,10,10,8,47,57,46,33,16,24,53,45,52,42,47,21,13,39,39,23,61,61,43,29,61,61,61,33,35,50,55,61,43,47,61,55,38,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,59,59,59,59,59,49,41,45,61,61,61,61,38,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,59,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,59,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,52,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,47,61,61,61,61,61,52,61,61,61,50,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,57,57,53,57,57,57,47,57,57,57,52,57,53,53,57,57,53,57,57,57,47,47,47,53,57,57,57,57,57,57,53,57,57,57,57,57,57,57,57,53,57,57,53,57,57,53,47,53,57,57,53,57,53,57,57,57,57,53,52,57,57,43,52,53,57,53,47,47,47,57,57,57,57,57,57,53,53,57,52,57,57,52,57,52,57,57,57,52,38,43,32,57,53,46,57,53,57,35,57,57,53,53,57,52,43,57,57,52,52,50,57,57,53,57,50,57,39,57,25,41,57,57,52,57,52,41,52,52,53,52,57,57,57,57,53,57,53,39,52,57,57,53,46,53,52,52,57,53,53,53,53,53,53,53,53,53,52,57,46,52,40,52,53,47,47,53,53,52,26,46,57,57,45,57,47,37,46,57,53,39,39,44,37,46,57,41,47,41,47,53,53,44,45,30,53,57,37,39,45,57,33,57,45,37,51,44,30,44,29,37,37,45,48,48,41,32,45,30,30,51,30,48,46,45,34,29,47,48,32,45,30,48,39,35,48,46,36,29,29,38,38,27,48,34,37,29,29,39,16,51,43,13,39,36,32,30,29,30,46,32,38,29,38,33,29,33,24,32,24,43,46,36,29,37,41,36,33,46,36,32,29,51,35,29,32,29,33,43,13,11]), \
		('HD-1-T7', \
			'NNNNNNNNNNNNNNNNNNTCTANNNTAATTTTGNTTAACTTTAAGAAGGAGATATACAGAATTCATGAAACATATGCCGAGAAAGATGTATAGTTGTGACTTTGAGACAACTACTAAAGTGGAAGACTGTAGGGTATGGGCGTATGGTTATATGAATATAGAAGATCACAGTGAGTACAAAATAGGTAATAGCCTGGATGAGTTTATGGCGTGGGTGTTGAAGGTACAAGCTGATCTATATTTCCATAACCTCAAATTTGACGGAGCTTTTATCACTAACTGGTTGGAACGTAATGGTTTTAAGTGGTCGGCTGACGGATTGCCAAACACATATAATACGATCATATCTCGCATGGGACAATGGTACATGATTGATATATGTTTAGGCTACAAAGGGAAACGTAAGATACATACAGTGATATATGACAGCTTAAAGAAACTACCGTTTCCTGTTAAGAAGATAGCTAAAGACTTTAAACTAACTGTTCTTAAAGGTGATATTGATTACCACAAAGAAAGACCAGTCGGCTATAAGATAACACCCGAAGAATACGCCTATATTAAAAACGATATTCAGATTATTGCGGAAGCTCTGTTAATTCAGTTTAAGCAAGGTTTAGACCGGATGACAGCAGGCAGTGACAGTCTAAAAGGTTTCAAGGATATTATAACCACTAAGAAATTCAAAAAGGTGTTTCCTACATTGAGTCTTGGACTCGATAAGGAAGTGAGATACGCCTATAGAGGTGGTTTTACATGTTAAATGATAGGTTCAAAGAAAAAGAAATCGGAGAAGGCATGGTCTTCGATGTTAATAGTCTATATCCTGCACAGATGTATAGCCGTCTCCTTCCATATGGTGAACCTATAGTATTCGAGGGTAAATACGTTTGGGACGAAGATTACCCACTACACATACAACATATCANATGTGAGTTCGAATTGAAAGAGGGCTATATACCCACTATACAGATAAAAGAAGTNNNTTTATAAAGGNAATGTGTACCTAAAAGTANCGGCGGGNNATANCGACNCTGGNNTCAANGTAGACNAGANNNGAAAGACACTACGATTANNTAACGTGANNNATCANCGGCTNAATTTNAAGCACTANNGGTTTNNNTTAANNNN', \
			[3,3,4,4,6,4,3,5,5,4,6,5,6,7,8,7,7,5,23,53,26,23,8,8,6,10,10,43,36,57,57,33,20,8,57,39,34,28,61,61,61,61,61,37,34,61,61,61,61,61,61,35,24,38,47,61,61,43,61,61,61,61,61,37,49,49,61,55,61,61,61,52,52,28,30,30,28,28,26,28,33,52,44,55,61,59,49,59,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,59,55,61,61,61,61,61,61,61,61,59,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,50,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,53,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,57,57,57,57,57,57,53,53,57,57,57,52,53,57,46,53,57,57,57,57,53,53,57,52,57,53,53,53,57,47,57,57,57,53,52,57,57,53,53,53,57,53,34,57,52,57,57,57,57,53,53,53,57,57,57,57,50,57,57,53,52,51,24,57,53,53,50,47,57,57,37,33,38,37,41,53,57,50,57,47,52,47,57,53,53,57,50,57,43,45,53,53,53,45,46,43,43,53,52,52,30,52,38,57,50,50,53,53,52,53,38,52,57,57,52,50,53,53,50,50,57,42,52,32,57,53,50,46,50,42,52,23,57,52,46,35,52,53,53,43,40,28,48,40,37,53,57,53,52,45,47,37,57,38,41,30,57,40,24,57,57,40,35,41,35,53,57,38,47,40,57,41,41,53,40,30,44,30,44,17,45,35,27,57,53,32,44,57,57,57,53,28,53,57,28,33,31,26,46,46,30,39,39,17,43,14,28,32,45,47,33,48,39,38,30,30,38,29,46,29,47,29,25,41,41,38,39,48,32,31,47,31,36,28,36,29,48,26,28,41,29,29,32,25,47,29,29,30,29,29,27,30,32,29,41,48,28,34,29,38,24,29,23,29,32,29,29,43,30,32,29,34,23,29,29,29,29,29,28,33,35,28,29,35,24,40,28,51,46,27,29,51,28,34,27,28,28,30,46,29,30,51,29,28,22,25,33,30,23,14,11,35,34,31,39,22,26,28,31,28,41,21,41,25,43,29,30,18,26,29,20,34,27,29,28,33,9,29,20,43,22,29,33,31,39,29,41,23,30,24,23,27,29,23,51,27,18,30,27,44,28,31,36,29,26,18,26,28,29,51,26,22,28,34,28,22,15,27,13,20,30,16,24,44,51,29,10,22,28,14,25,5,6,7,51,51,15,16,28,15,36,12,29,28,8,14,13,17,29,14,30,23,22,16,20,32,32,44,39,34,13,21,27,9,30,23,19,20,24,44,19,8,5,13,12,16,9,28,11,22,31,9,21,22,31,24,6,4,14,21,19,19,9,18,10,25,12,18,29,7,12,12,16,6,5,7,26,10,36,20,20,30,12,10,12,15,21,29,18,31,30,25,10,6,8,16,11,24,25,12,23,10,10,4,4,6,10,24,24,13,5,16,13,18,16,20,7,31,17,12,26,14,4,25,15,16,11,14,12,13,11,9,3,20,19,12,31,16,5,4,4,30,16,13,23,5,6,4,5]), \
		('HD-1-T7-Term','NNNNNNNNNNNNNNNTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGGTGATGATGATGACCGCTACCGCCGTTGATTGTGAATGTGTCATCAACCAGAACCACCCCGCCCGGCACTTGCGCAGGTTTAGGCTTCATTTTCCGACTGAATCCGACTTTGAAATTCTCAAACGTAACCTCTTTCTTAATCTTGTCGGTCATTCCCGCACATTTAACACTAAATTTTATATCAGTGTAATCATCTGGACTACCTTCTACTAACTTACCATCTACTTCTTTCATATAGATGTCTTGTATATAGGTCTTCTGTCTCAGATATTTAGCTCTTTTGAATGTACTTTCATGTGCCCAGTATCCCAATTTCTTAGGGTCAACTATATCTTTTATTACATCAGGTATCTCTGTACCCGTTAAATGTATGCTGTCAGTATCACATTATATTATCCGATCATAACAAGCCTGTGCCGCTGTAATTGTCGTGTATCTAGCCCATGCAGTGATGAAAACGCCCATAGGCGTATAAACAGGGTCTTTTGTTTCCTCTTCTCCAAGTCTGAAACCTAGCGCCCCATTCTCTTTTAAATAAAGGACTTTCCCTGTAACATCAGGGTTACTAGCGAATTTACCGTATAGACTGTTTAACATCAGTTTTGCTAGTTGCTTGATCGCTCCTACTGATGTCGTCATGATGTACGTCCATTTATCTATAAAATCTTTAAACAAACCTGTAGTTGCTTTAAATTTTAAGCCGCTGATATATTCAACGTTATATAAATCGTAGTGTTCTTTCATTAATTCTANGTCTACATTTGACAACCAGAGGTCGGCTATCTCCCCGCCGCTACTTTTTAGGTACACATTN', \
			[1,1,2,1,6,4,5,3,4,7,5,7,6,7,7,20,46,29,29,22,25,57,30,51,41,31,37,52,21,29,23,24,21,27,33,57,57,57,19,23,53,57,31,21,61,61,61,42,47,61,43,40,61,47,47,61,47,47,61,47,61,61,47,47,61,47,38,61,49,49,61,39,36,43,36,23,19,18,12,16,13,44,44,55,52,61,61,61,61,61,59,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,61,61,61,52,55,55,43,43,50,37,50,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,42,61,61,61,61,61,61,61,61,61,61,61,42,61,61,61,61,53,61,52,61,47,61,61,61,53,57,57,57,57,57,53,57,57,57,47,57,52,53,57,57,57,53,53,57,53,57,53,53,57,57,57,50,57,57,57,53,57,57,53,57,53,52,57,53,53,47,42,37,47,57,53,53,53,47,53,57,57,50,45,57,57,57,57,53,47,57,53,33,57,53,57,57,39,57,47,47,41,47,53,53,53,53,53,53,57,53,57,57,38,47,47,53,52,57,52,57,53,52,52,52,57,57,57,43,41,52,52,52,57,57,53,45,57,57,57,53,44,52,52,53,53,53,57,47,53,50,53,57,41,38,57,39,53,37,36,53,42,52,57,53,32,57,50,53,38,52,53,37,47,36,53,57,47,53,40,57,47,53,53,52,46,41,38,41,53,53,57,38,34,53,53,33,37,52,44,30,38,39,53,40,53,32,52,52,52,22,46,57,53,34,57,40,47,46,34,32,30,34,33,39,44,38,37,45,40,40,24,47,37,48,48,53,35,28,44,37,41,30,31,47,31,47,44,30,21,32,47,30,30,27,46,22,44,51,20,30,44,29,46,40,32,30,32,30,29,38,39,46,45,29,44,9,46,29,29,29,29,29,24,36,51,41,30,29,26,28,29,24,43,29,22,43,13,51,29,29,36,29,30,43,41,41,41,36,27,39,26,39,29,31,32,29,30,29,29,29,36,51,39,39,29,29,15,32,27,29,30,24,29,21,26,12,9]), \
		('HD-2-F0_303_F','NNNNNNNNNNNTNNGCTNNNGGGAACGTAAGANACATACAGTGATATAAGACAGCTTAAAGAAACTACCGTTTCCTGTTAAGAAGATAGCTAAAGACTTTAAACTAACTGTTCTTAAAGGTGATATTGATTACCACAAAGAAAGACCAGTCGGCTATAAGATAATACCCGAAGAATACGCCTATATTAGAAACGATATTCAGATTATTGCGGAAGCTCTGTTAATTCAGTTTAAGCAAGGTTTAGACCGGATGACAGCAGGCAGTGACAGTCTAAAAGGTTTCAAGGATATTATAACCACTAAGAAATTCAAAAAGGTGTTTCCTACATTGAGTCCTGGACTCGATAAGGAAGTAAGATACGCCTATAAAGGTGGTTTTACATGGTTAAATGATAGGTTCAAAGAAAAAGAAATCGGAGAGGGCATAGTCTTCGATGTTAATAGTCTATATCCTGCACAGATGTATAGCCGCCTCCTTCCATATGGTGAACCTATAGTATTCGAGGGTAAATACGTTTGGGACGAAGATTACCCACTACACATACAGCATATCAGATGTGAGTTCGAATTGAAAGAGGGCTATATACCCAATATACAGATAAAAAGAAGTAGGTTTTATAAAGGTAATGAGTACCTAAAAAGTAGCGGCGGGGAGATAGCCGACCTCTGGTTGTCAAATGTAGACCTAGAATTAATGAAAGAACACTACGATTTATATAACGTTGAATATATCAGCGGCTTAAAATTTAAAGCAACTACAGGTTTGTTTAAAGATTTTATAGATAAATGGACGTACATCAAGACGACATCAGAAGGAGCGATCAAGCAACTAGCAAAACTGATGTTAAACAGTCTATACGGTAAATTCGCTAGTAACCCTGATGTTACAGTGAAAGTCCCTTATTTAAAAGAGAATGGGGCGCTAGGTTTCAGACTTGGAGAAGAGGAAACAAAAGACCCTGTTTATACNCCTATGGGCGTTTTCATCANNN', \
			[1,1,3,3,4,6,5,7,5,5,6,10,7,7,25,15,13,8,7,8,10,43,41,33,11,30,53,50,47,42,53,20,7,30,42,22,61,61,43,24,61,61,61,33,40,43,55,55,47,49,61,52,49,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,61,61,61,59,59,59,55,52,59,49,41,36,59,59,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,44,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,59,61,61,61,61,61,61,61,61,61,61,61,55,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,53,61,61,61,61,61,61,61,61,61,50,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,57,57,57,57,57,53,57,57,57,57,57,57,57,57,57,53,52,57,57,57,57,57,57,47,37,47,53,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,53,53,57,53,42,40,57,47,47,57,57,53,57,57,57,57,52,57,57,52,57,57,53,57,53,41,47,57,57,57,57,57,57,52,52,57,50,57,57,57,50,57,57,57,57,52,32,43,29,57,53,41,43,53,52,27,33,41,53,53,43,41,36,53,57,34,53,38,57,57,53,57,42,57,50,52,47,57,53,53,50,57,40,52,53,39,52,40,57,53,57,57,43,57,37,53,50,52,57,53,57,53,57,42,57,53,53,53,53,53,40,57,53,53,39,57,46,41,44,37,53,47,53,53,53,40,47,57,57,46,37,57,53,37,41,57,37,53,57,52,34,38,57,46,40,33,35,40,43,41,53,30,53,57,31,39,37,57,38,57,53,37,57,48,37,47,27,51,34,40,47,47,41,37,35,40,47,51,41,47,46,44,32,29,38,47,30,44,21,46,32,32,37,40,44,24,29,46,30,27,47,33,32,30,51,48,27,51,37,20,37,29,29,29,28,32,48,29,47,29,46,44,36,30,25,36,29,51,51,32,29,44,47,31,28,46,28,36,29,51,30,44,39,36,29,43,29,36,29,29,41,18,43,23,43,51,30,29,43,28,28,29,29,35,34,29,29,35,29,39,36,29,41,35,29,43,17,29,23,29,28,29,29,39,39,51,32,28,29,27,51,41,30,36,29,26,51,36,29,51,51,39,22,46,25,31,30,28,41,46,46,39,27,29,25,29,10,12,16,16,44,29,28,24,14,36,17,22,29,33,34,31,16,32,33,12,46,22,34,31,51,28,21,26,46,51,39,22,43,18,36,20,22,38,22,51,26,22,18,23,13,7,11,17,26,26,28,31,44,38,29,31,28,51,44,17,28,27,29,28,28,6,5,5]), \
		('HD-2-T7', \
			'NNNNNNNNNNNNNNTNNNNCTAGAATAATTTTGNTTAACTTTAAGAAGGAGATATACAGAATTCATGAAACATATGCCGAGAAAGATGTATAGTTGTGACTTTGAGACAACTACTAAAGTGGAAGACTGTAGGGTACGGGCGTATGGTTATATGAATATAGAAGATCACAGTGAGTACAAAATAGGTAATAGCCTGGATGAGTTTATGGCGTGGGTGTTGAAGGTACAAGCTGATCTATATTTCCATAACCTCAAATTTGACGGAGCTTTTATCATTAACTGGTTGGAACGTAATGGTTTTAAGTGGTCGGCTGACGGATTGCCAAACACATATAATACGATCATATCTCGCATGGGACAATGGTACATGATTGATATATGTTTAGGCTACAAAGGGAAACGTAAGATACATACAGTGATATAAGACAGCTTAAAGAAACTACCGTTTCCTGTTAAGAAGATAGCTAAAGACTTTAAACTAACTGTTCTTAAAGGTGATATTGATTACCACAAAGAAAGACCAGTCGGCTATAAGATAATACCCGAAGAATACGCCTATATTAGAAACGATATTCAGATTATTGCGGAAGCTCTGTTAATTCAGTTTAAGCAAGGTTTAGACCGGATGACAGCAGGCAGTGACAGTCTAAAAGGTTTCAAGGATATTATAACCACTAAGAAATTCAAAAAGGTGTTTCCTACATTGAGTCCTGGACTCGATAAGGAAGTAAGATACGCCTATAAAGGTGGTTTTACATGGTTAAATGATAGGTTCAAAGAAAAAGAAATCGGAGAGGGCATAGTCTTCGATGTTAATAGTCTATATCCTGCACAGATGTATAGCCGCCTCCTTCCATATGNN', \
			[2,1,1,3,3,4,4,5,5,4,6,4,5,7,10,7,7,6,7,52,19,23,24,43,39,24,30,26,43,57,57,33,21,9,57,24,47,28,61,61,50,29,50,43,47,61,47,61,61,61,61,35,24,39,47,61,61,43,61,61,61,61,61,47,61,49,61,61,61,61,61,59,41,27,23,23,27,31,20,27,20,36,55,59,59,59,49,59,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,59,61,59,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,44,59,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,61,61,53,61,53,61,61,61,61,61,61,61,61,61,61,61,50,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,52,61,61,61,61,61,61,61,61,61,61,53,61,53,61,61,61,61,61,61,61,61,43,61,61,61,61,61,61,61,61,61,61,61,53,61,61,61,40,61,61,61,57,53,57,57,57,57,57,53,57,57,53,57,53,53,57,57,57,57,57,57,53,57,57,53,57,57,57,57,53,47,47,53,53,57,57,47,57,57,38,53,53,47,57,57,53,57,57,57,57,53,57,57,53,57,57,53,47,46,57,52,57,53,45,52,52,57,53,53,57,52,57,47,47,52,57,45,57,53,57,45,57,52,57,57,57,41,48,57,52,43,50,52,57,53,31,57,50,53,45,53,57,41,37,50,53,46,57,53,47,41,41,52,41,43,35,57,53,53,52,43,53,57,37,36,52,47,53,53,53,52,53,46,38,35,53,53,53,53,42,57,57,33,40,57,37,35,57,57,57,46,39,46,27,57,44,57,53,44,57,41,52,45,30,46,38,46,53,40,53,42,53,44,57,57,43,27,57,41,38,53,41,31,53,24,43,28,57,33,37,53,53,36,40,28,41,34,35,41,53,47,53,34,23,51,45,39,45,28,46,45,30,46,46,45,30,31,29,48,51,36,25,46,31,51,41,33,46,32,29,47,28,51,19,30,29,29,46,43,29,46,51,51,46,51,28,46,51,35,41,33,29,29,46,22,38,16,46,45,25,29,43,29,27,27,43,33,27,26,26,31,28,43,30,29,29,36,32,35,34,35,36,30,29,32,33,36,29,43,38,43,25,25,36,28,28,36,20,40,29,29,31,35,34,29,39,23,23,32,34,20,23,20,20,13,29,23,29,29,29,27,8,5]), \
		('HD-2-T7-Term', \
			'NNNNNNNNNCNNCTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGCTCGAGTTAGTGGTGGTGATGATGATGACCGCTACCGCCTTTGATTGTGAATGTGTCATCAACCAGAACCACCCCGCCCGGCACTTGCACAGGCTTAGGCTTCATTTTCCGACTGAATCCGACTTTGAAATTCTCAAACGTAACCTCTTTCTTAATCTTGTCAGTCATTCCCGCACATTTAACACTAAATTTTATATCAGTGTAATCATCTGGACTACCTTCTACTAACTTACCATCTACTTCTTTCATATAGATGTCTTGTATATAGGTCTTCTGTCTCAGATATTCAGCTCTTTTGAATGTACTTTCATGTGCCCAGTATCCCAATTTCTTAGGGTCAACTATATCTTTTATTACATCTGGTATCTCTGTACCCGTTAAATGTATGCTGTCAGTATCACAGTATATTATCCGATCATAACAAGCCTGTGCCGCTGTAATTGTCGTGTATCTAGCCCATGCAGTGATGAAAACGCCCATAGGTGTATAAACAGGGTCTTTTGTTTCCTCTTCTCCAAGTCTGAAACCTAGCGCCCCATTCTCTTTTAAATAAGGGACTTTCACTGTAACATCAGGGTTACTAGCGAATTTACCGTATAGACTGTTTAACATCAGTTTTGCTAGTTGCTTGATCGCTCCTTCTGATGTCGTCTTGATGTACGTCCATTTATCTATAAAATCTTTAAACAAACCTGTAGTTGCTTTAAATTTTAAGCCGCTGATATATTCAACGTTATATAAATCGTAGTGTTCTTTCATTAATTCTAGGTCTACATTTGACAACCAGAGGTCGGCTATCTCCCCGCCGCTACTTTTTAGGTACTCATTACCTTTATAAAACCTACTTCTTTTTATCTGTATATTGGGTATATAGCCCTCN', \
			[4,1,1,4,3,3,4,5,4,10,6,8,12,10,33,18,26,32,34,53,24,57,57,27,38,52,19,29,22,29,16,17,40,57,57,57,20,24,61,61,33,24,61,61,61,43,47,61,47,41,61,47,47,61,47,55,61,61,61,61,61,61,61,61,61,61,61,55,49,36,43,59,37,27,59,41,41,43,59,59,61,61,59,61,61,59,61,61,61,61,61,61,59,61,59,61,61,59,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,59,44,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,55,44,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,55,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,59,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,53,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,50,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,47,61,61,61,61,61,61,61,57,53,57,47,57,57,47,57,53,57,47,47,57,53,57,57,57,57,57,57,57,57,53,57,57,57,57,53,57,47,57,53,57,50,57,57,57,57,52,57,57,57,53,57,57,53,53,57,57,57,53,57,57,57,47,47,47,57,57,53,57,53,57,52,45,57,53,52,41,33,57,53,57,57,57,53,57,53,57,57,45,57,52,53,53,53,53,57,57,53,50,57,57,57,57,57,50,53,52,50,53,57,57,53,57,53,57,52,47,53,52,52,57,53,57,53,53,50,57,52,53,30,57,39,53,50,50,57,57,57,57,52,40,52,52,40,52,50,53,53,46,38,53,57,46,53,52,53,47,53,53,52,57,41,44,42,53,57,57,46,38,33,38,46,40,52,44,41,53,53,57,40,53,34,52,52,52,41,57,57,41,40,57,41,53,57,57,40,33,53,30,34,57,40,41,53,45,43,29,33,45,39,31,40,41,25,44,53,47,47,45,44,47,44,48,47,42,37,45,44,41,30,51,29,30,46,30,40,43,48,51,40,31,30,40,32,31,44,30,39,29,29,46,39,48,29,30,41,29,29,29,43,51,43,36,31,30,29,31,32,37,29,37,38,13,46,29,35,29,29,29,29,29,29,29,36,29,41,30,26,29,32,37,29,29,29,28,29,36,46,51,46,41,23,24,15,18,30,28,38,29,25,30,31,29,21,32,30,51,28,33,43,31,51,51,26,13,26,28,29,29,27,41,35,31,51,46,51,28,29,29,28,29,46,29,29,28,31,24,33,36,51,34,28,29,29,27,29,40,28,33,36,21,33,20,8])]

	def tearDown(self):
		self.a = 1
		# self.widget.dispose()
