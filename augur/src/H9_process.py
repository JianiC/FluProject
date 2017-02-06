import time, argparse,re,os, socket
import matplotlib as mpl
if socket.gethostname() not in ['olt', 'rneher-iMac']:
	mpl.use('pdf')
from virus_filter import flu_filter, fix_name
from virus_clean import virus_clean
from tree_refine import tree_refine
from tree_titer import HI_tree
from fitness_model import fitness_model
from process import process, virus_config
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import numpy as np
from itertools import izip

virus_config.update({
	# data source and sequence parsing/cleaning/processing
	'virus':'H9',
	'alignment_file':'/Users/yujiazhou/Documents/nextflu/H9_nextflu-master/augur/src/data/H9_gisaid_epiflu_sequence.fasta',
	'outgroup':'A/duck/HongKong/147/1977',
	#'force_include':'/Users/yujiazhou/Documents/nextflu/H9_nextflu-master/augur/src/data/H9_HI_strains.txt',
	'force_include_all':False,
	'date_spec':'year',
	'max_global':True,   # sample as evenly as possible from different geographic regions
	#'max_globalh':True,
	'cds':[0,None], # define the HA1 start i n 0 numbering
	'n_iqd':5,
	'min_mutation_frequency':0.01,
	# define relevant clades in canonical HA1 numbering (+1)
	# numbering starting at HA1 start, adding sp to obtain numbering from methionine
	'''
	'clade_designations': { "Y439":[('HA1',122,'F'), ('HA1',353,'P')],
							"Korea":[('HA1',107,'M'), ('HA1',122,'F'), ('HA1',127,'R'), ('HA1',130,'K'), ('HA1',132,'L'), ('HA1',134,'L'), ('HA1',179,'D'), ('HA1',212,'I'), ('HA1',299,'T'), ('HA1',353,'P'), ('HA1',473,'K')],
						   	"G1":[('HA1',353,'P'), ('HA1',473,'K')],
							"Ck-Bei":[('HA1',107,'M'), ('HA1',299,'T'), ('HA1',473,'K')],
						   	"G9":[('HA1',107,'M'), ('HA1',299,'T'), ('HA1',473,'K')],
						   	"Y280":[('HA1',299,'T'), ('HA1',473,'K')]
							},
	'''
	#'epitope_masks_fname':'/Users/yujiazhou/Documents/nextflu/H9_nextflu-master/augur/source-data/H9_epitope_masks.tsv',
	#'epitope_mask_version':'wolf',
	#'HI_fname':'/Users/yujiazhou/Documents/nextflu/H9_nextflu-master/augur/src/data/H9_HI_titers.txt',
	'auspice_prefix':'H9_',
	'html_vars': {'coloring': 'ep, ne, rb, lbi, dfreq, region, date, cHI',
				   'gtplaceholder': 'HA1 positions...',
					'freqdefault': 'Y439, Korea, G1, Ck-Bei, G9, Y280'},
	'js_vars': {'LBItau': 0.0005, 'LBItime_window': 0.5, 'dfreq_dn':2},
	'excluded_tables': ['NIMR_Sep2012_08.csv'], #, 'nimr-sep-2010-table8', 'nimr-sep-2010-table8','NIMR_Sep2012_11.csv'],
	'layout':'auspice',
	'min_aamuts': 1,
#	'predictors': ['dfreq', 'cHI']												# estimate
	'predictors': { 'dfreq': [2.50, 2.84], 'cHI': [1.68, 0.45] }				# fix predictor: [value, std deviation]
	})


class H9_filter(flu_filter):
	def __init__(self,min_length = 0, **kwargs):
		'''
		parameters
		min_length  -- minimal length for a sequence to be acceptable
		'''
		flu_filter.__init__(self, **kwargs)
		self.min_length = min_length
		self.vaccine_strains =[
				{
					"strain": "A/Chicken/HongKong/G9/97",
					"db": "GISAID",
					"accession": "EPI_ISL_1263",
					"date": "1997-07-01",
					"seq": "ATGGAAATAATAGCACTAATAGCTATACTGGTAGTGACAAAAACAAGCAATGCAGATAAAATTTGCATTGGCTACCAGTCAACAAACTCCACAGAAACTGTTGATACACTAGTAGAAAACAATGTCCCTGTGACACATACCAAAGAATTGCTCCACACAGAGCACAATGGAATGCTATGTGCAACAAACCTGGGGCACCCTCTCATCCTAGACACCTGCACCATCGAAGGGTTGGTGTACGGCAACCCTTCCTGTGATTTGCTACTGGGAGGGAAAGAATGGTCTTACATTGTCGAAAGATCATCAGCTGTCAATGGGATGTGTTACCCTGGAAGGGTAGAGAACCTGGAAGAACTCAGGTCTTTTTTCAGCTCCGCTCGCTCCTACAAAAGACTCCTGCTCTTTCCAGACAGAACTTGGAATGTGACTTACACTGGGACAAGCAAAGCATGTTCAAACTCATTCTACAGAAGTATGAGATGGCTGACACACAAGAGCGATTCTTACCCTATTCAAGACGCCCAATATACTAACGATTGGGGAAAGAATATTCTCTTCATGTGGGGCATACACCACCCACCTACTGATACTGAGCAAATAAATCTATACAAAAAAGCTGATACAACAACAAGTATAACAACGGAAGATATCAATCGAACTTTCAAACCAGTGATAGGGCCAAGGCCTCTTGTCAATGGTCAACAAGGGAGAATTGATTATTATTGGTCAGTACTAAAGCCAGGCCAGACACTGCGAGTGAGATCCAATGGGAATCTAATTGCCCCATGGTATGGACACATTCTTTCAGGAGAAAGCCATGGAAGAATCTTGAAGACCGATTTGAGTAGTGGCAACTGCGTAGTACAATGCCAAACTGAGAAAGGTGGTTTGAACACGACCTTGCCATTCCACAATGTCAGCAAGTATGCATTTGGGAACTGCCCCAAATATGTTGGAGTGAAGAGTCTCAAACTGGCAGTTGGTCTAAGGAATGTTCCTGCTGCATCATATAGAGGGCTCTTCGGTGCCATAGCTGGATTCATAGAAGGCGGTTGGCCAGGACTAGTTGCAGGCTGGTACGGGTTTCAGCATTCAAATGATCAAGGGGTTGGAATGGCCGCAGATAGGGAATCAACTCAAGAAGCAGTTGACAAGATAACATCCAAAGTAAATAACATAATCGACAAAATGAACAAGCAGTATGGA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T--------------------------------------------------------------------------------------"}
			]
		tmp_outgroup = SeqIO.read('/Users/yujiazhou/Documents/nextflu/H9_nextflu-master/augur/source-data/H9_outgroup.gb', 'genbank')
		genome_annotation = tmp_outgroup.features
		self.cds = {x.qualifiers['gene'][0]:x for x in genome_annotation
				if 'gene' in x.qualifiers and x.type=='CDS' and
				x.qualifiers['gene'][0] in ['SigPep', 'HA1', 'HA2']}
		self.outgroup = {
			'strain': 'A/duck/HongKong/147/1977',
			'db': 'IRD',
			'accession': 'AY206671',
			'date': '2003-03-03',
			'country': 'HongKong',
			'region': 'EastAsia',
			'seq': str(tmp_outgroup.seq).upper()
		}

class H9_clean(virus_clean):
	def __init__(self,**kwargs):
		virus_clean.__init__(self, **kwargs)

	def clean_outbreaks(self):
		"""Remove duplicate strains, where the geographic location, date of sampling and sequence are identical"""
		virus_hashes = set()
		new_viruses = []
		for v in self.viruses:
			try:
				geo = re.search(r'A/([^/]+)/([^/]+)/', v.strain).group(2)
			except:
				print "clean outbreaks:, couldn't parse geo of ",v.strain
				continue
			if geo:
				vhash = (geo, v.date, str(v.seq))
				if vhash not in virus_hashes:
					new_viruses.append(v)
					virus_hashes.add(vhash)

		self.viruses = MultipleSeqAlignment(new_viruses)

	'''def clean_reassortants(self):
		from seq_util import hamming_distance as distance
		#Remove viruses from the outbreak of triple reassortant pH1N1
		remove_viruses = []

		reassortant_seqs = [
			"ATGAAGACTATCATTGCTTTTAGCTGCATTTTATGTCTGATTTTCGCTCAAAAACTTCCCGGAAGTGACAACAGCATGGCAACGCTGTGCCTGGGACACCATGCAGTGCCAAACGGAACATTAGTGAAAACAATCACGGATGACCAAATTGAAGTGACTAATGCTACTGAGCTGGTCCAGAGTTCCTCAACAGGTGGAATATGCAACAGTCCTCACCAAATCCTTGATGGGAAAAATTGCACACTGATAGATGCTCTATTGGGGGACCCTCATTGTGATGACTTCCAAAACAAGGAATGGGACCTTTTTGTTGAACGAAGCACAGCCTACAGCAACTGTTACCCTTATTACGTGCCGGATTATGCCACCCTTAGATCATTAGTTGCCTCATCCGGCAACCTGGAATTTACCCAAGAAAGCTTCAATTGGACTGGAGTTGCTCAAGGCGGATCAAGCTATGCCTGCAGAAGGGGATCTGTTAACAGTTTCTTTAGTAGATTGAATTGGTTGTATAACTTGAATTACAAGTATCCAGAGCAGAACGTAACTATGCCAAACAATGACAAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACGGACAAGGACCAAACCAACCTATATGTCCAAGCATCAGGGAGAGTTATAGTCTCTACCAAAAGAAGCCAACAAACTGTAATCCCGAATATCGGGTCTAGACCCTGGGTAAGGGGTGTCTCCAGCATAATAAGCATCTATTGGACGATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCCCCTCGGGGTTACTTCAAAATACAAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACACATTGATGAATGCAATTCTGAATGCATTACTCCAAATGGAAGCATTCCCAATGACAAACCTTTTCAAAATGTAAACAAGATCACATATGGAGCCTGTCCCAGATATGTTAAGCAAAACACCCTGAAATTGGCAACAGGAATGCGGAATGTACCAGAGAAACAAACTAGAGGCATATTCGGCGCAATTGCAGGTTTCATAGAAAATGGTTGGGAGGGAATGGTAGACGGTTGGTACGGTTTCAGGCATCAGAATTCTGAAGGCACAGGACAAGCAGCAGATCTTAAAAGCACTCAAGCAGCAATCAACCAAATCACCGGGAAACTAAATAGAGTAATCAAGAAAACAAACGAGAAATTCCATCAAATCGAAAAAGAATTCTCAGAAGTAGAAGGAAGAATTCAGGACCTAGAGAAATACGTTGAAGACACTAAAATAGATCTCTGGTCTTACAACGCTGAGATTCTTGTTGCCCTGGAGAACCAACATACAATTGATTTAACCGACTCAGAGATGAGCAAACTGTTCGAAAGAACAAGAAGGCAACTGCGGGAAAATGCTGAGGACATGGGCAATGGTTGCTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATCAGAAATGGAACTTATGACCATGATATATACAGAAACGAGGCATTAAACAATCGGTTCCAGATCAAAGGTGTTCAGCTAAAGTCAGGATACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGCTTTTTGCTTTGTGTTGTTCTGCTGGGGTTCATTATGTGGGCCTGCCAAAAAGGCAACATTAGGTGCAACATTTGCATTTGA",
			"ATGAAGACTATCATTGCTTTTAGCTGCATCTTATGTCAGATCTCCGCTCAAAAACTCCCCGGAAGTGACAACAGCATGGCAACGCTGTGCCTGGGGCATCACGCAGTACCAAACGGAACGTTAGTGAAAACAATAACAGATGACCAAATTGAAGTGACTAATGCTACTGAGCTGGTCCAGAGTACCTCAAAAGGTGAAATATGCAGTAGTCCTCACCAAATCCTTGATGGAAAAAATTGTACACTGATAGATGCTCTATTGGGAGACCCTCATTGTGATGACTTCCAAAACAAGAAATGGGACCTTTTTGTTGAACGAAGCACAGCTTACAGCAACTGTTACCCTTATTATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACCCTGGAATTTACTCAAGAAAGCTTCAATTGGACTGGGGTTGCTCAAGACGGAGCAAGCTATTCTTGCAGAAGGGAATCTGAAAACAGTTTCTTTAGTAGATTGAATTGGTTATATAGTTTGAATTACAAATATCCAGCGCTGAACGTAACTATGCCAAACAATGACAAATTTGACAAATTGTACATTTGGGGGGTACACCACCCGGGTACGGACAAGGACCAAACCAGTCTATATATTCAAGCATCAGGGAGAGTTACAGTCTCCACCAAATGGAGCCAACAAACTGTAATCCCGAATATCGGGTCTAGACCCTGGATAAGGGGTGTCTCCAGCATAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCCCCTCGGGGTTACTTCAAAATACAAAGTGGGAAAAGCTCAATAATGAGGTCAGATGCACACATTGGCAACTGCAACTCTGAATGCATTACCCCAAATGGAAGCATTCCCAACGACAAACCTTTTCAAAATGTAAACAGAATAACATATGGGGCCTGTCCCAGATATGTTAAGCAAAACACTCTGAAATTAGCAACAGGAATGCGGAATGTACCAGAGAAACAAACTAGAGGCATATTCGGCGCAATCGCAGGTTTCATAGAAAATGGTTGGGAAGGGATGGTGGACGGTTGGTATGGTTTCAGGCATCAAAACTCTGAAGGCACAGGGCAAGCAGCAGATCTTAAAAGCACTCAAGCGGCAATCAACCAAATCACCGGGAAACTAAATAGAGTAATCAAGAAGACGAATGAAAAATTCCATCAGATCGAAAAAGAATTCTCAGAAGTAGAAGGGAGAATTCAGGACCTAGAGAGATACGTTGAAGACACTAAAATAGACCTCTGGTCTTACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATTTAACTGACTCAGAAATGAACAAACTGTTCGAAAGGACAAGGAAGCAACTGCGGGAAAATGCTGAGGACATGGGCAATGGATGCTTTAAAATATATCACAAATGTGACAATGCCTGCATAGGATCAATCAGAAATGGAACTTATGACCATGATGTATACAGAGACGAAGCAGTAAACAATCGGTTCCAGATCAAAGGTGTTCAGCTGAAGTTAGGATACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGCTTTTTGCTTTGTGCTGTTCTGCTAGGATTCATTATGTGGGCATGCCAAAAAGGCAACATTAGGTGCAACATTTGCATTTGA",
			"ATGAAGACTAGTAGTTCTGCTATATACATTGCAA------------------------CCGCAAATG---------CAGACACATTATGTATAGGTTATCATGCAGTACTAGAAAAGAATGTAACAGTAACACACTCTGTTAACCAAACTGAGAGGGGTAGCCCCATTGCATTTG--------------------GGTAAATGTAACATTGCTGGCTGGATCC------------------------------------TGGGAAATCCAGAGTGTGACACTCTCCACAGCAAGCTCATGGTCCTACATCGTGGAAACATCTAAGACAATGGAACGTGCTACCCAGGAGATTTCATCAATTATGAGGAGCTAAGGTCATCATTTGAAAGGTTTGAGATATTACAAGTTCATGGCCCAATCATGACTCGAACAAAGGTTCCTCAAGCTGGAGCAA---------------------------AAAGCTTCTACAAAAATTTAATATGGCTAGTTAAAAAAGGAAATTCATACCCAA------------------------------AGCTCAGCAAATCCTACATTTGGGGCATTCACCATCCATCTACTAGTGCTGACCAA-------CAAAGTCTCTATCAGAGTGCAGATGCATATGTTTTATCAAAATACAGCAAGAAGTTCAAG--CCGGAAATAGCAGTAAGACCCAAAGTGAGGGATCAAGAAGGGAGAATGAACTATTACTGGACACTAGTAGAGCCGGGAGACAAAATAACATTCGAAGCAACTGGAAATCTATTGGTACCGAGATATGCATTCGCAATGGAAA----GAAATGCTGGATTATCATTTCAGATACACCAGTCCACGATTGCAATACAACTTGTCAGACACCCAAGGGTGCTATAAACACCAGCCTCCCATTTCAGAATATACATCCGATCACAATTGGAAAATGTCCCAAATATGTAAAAAGCACAAAATTGAGACTGGCCACAGGATTGAGGAATGTCCCGTCTATTCAATCTAGAGGCCTATTTGGGGCCATTGCCGGTTTCATTGAAGGGGGGTGGACAGGGATGGTAGATGGATGGTACGGTTATCACCATCAAAATGCGCAGGGGTCAGGATATGCAGCCGACCTGAAGAGCACACAGAATGCCATTGACAAGATTACTAACAAAGTAAATTCTGTTATTGAAAAGATGAATACACAGTTCACAGCAGTAGGTAAAGAGTTCAACCACCTGGAAAAAAGAATAGAGAATTTAAATAAAAAAGTTGATGATGGTTTCCTGGACATTTGGACTTACAATGCCGAACTGTTGGTTCTATTGGAAAATGAAAGAACTTTGGACTACCACGATTCAAATGTGAAAAACTTATATGAAAAGGTAAGAAGCCAGTTAAAAAACAATGCCAAGGAAATTGGAAACGGCTGCTTTGAATTTTACCACAAATGCGATAACACGTGCATGGAAAGTGTCAAAAATGGGACTTATGACTACCCAAAATACTCAGAGGAAGCAAAATTAAACAGAGAAGAAATAGATGGGGTAAAGCTGGAATCAACAAGGATTTACCAGATTTTGGCGATCTATTCAACTGTCGCCAGTTCATTGGTACTGGTAGTCTCCCTGGGGGCAATCATCTGGATGTGCTCTAATGGGTCTCTACAGTGTAGAATATGTATTTAA",
			"ATGAAGACTATCATTGCTTTGAGCTACATTTTATGTCTGGTTTTCGCTCAAAAACTTCCCGGAAATGACAACAGCACGGCAACGCTGTGCCTGGGGCACCATGCAGTGCCAAACGGAACGCTAGTGAAAACAATCACGAATGACCAAATTGAAGTAACTAATGCTACTGAGCTGGTTCAGAGTTCCTCAACAGGTAGAATATGCGACAGTCCTCACCAAATCCTTGATGGAGAAAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCATTGTGATGGCTTCCAAAACAAGGAATGGGACCTTTTTGTTGAACGCAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGTCTCCCTTAGGTCACTAGTTGCCTCATCAGGCACGCTGGAGTTTAACAATGAAAGCTTCAATTGGACTGGAGTCGCTCAGAATGGAACAAGCTCTGCTTGCAAAAGGAGATCCGATAAAAGTTTCTTTAGTAGATTGAATTGGTTGCACCAATTAAAATACAAATATCCAGCACTGAACGTGACTATGCCAAACAATGAAAAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACAGACAGTGACCAAATCAGCCTATATGCTCAAGCATCAGGGAGAGTCACAGTCTCTACCAAAAGAAGCCAACAAACTGTAATCCCGAATATCGGATCTGGACCCTGGGTAAGGGGTGTCTCCAGCAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTCGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGGTCAGATGCACCCATTGGCAAATGCAATTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGGCAAACCATTTCAAAATGTAAACAGGATCACATATGGGGCCTGTCCCAGATATGTTAAGCAAAACACTCTGAAATTGGCAACAGGGATGCGGAATGTGCCAGAGAAACAAACTAGAGGCATATTCGGTGCAATCGCGGGCTTCATAGAAAATGGTTGGGAGGGAATGATGGACGGTTGGTACGGTTTCAGGCATCAGAATTCTGAGGGCACAGGGCAAGCAGCAGATCTTAAAAGCACTCAAGCAGCAATCAACCAAATCAACGGGAAACTGAATAGGTTAATCGAGAAAACGAACGAGAAATTCCATCAAATTGAAAAAGAATTCTCAGAAGTAGAAGGGAGAATTCAGGACCTCGAGAAATATGTCGAGGACACTAAAATAGATCTCTGGTCGTACAATGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATCTAACTGACTCAGAAATGAACAAACTGTTTGAAAGAACAAAGAAGCAACTGAGGGAAAATGCTGAGGATATGGGCAATGGTTGTTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGGTCAATCAGAAATGGAACTTATGACCATGATGTATACAGAGACGAAGCATTGAACAACCGGTTCCAGATCAAAGGTGTTGAGCTGAAGTCAGGATACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGTTTTTTGCTTTGTATTGTTTTACTGGGGTTCATCATGTGGGCCTGCCAAAAAGGCAACATTAGGTGCAACATTTGCATTTGA",
			"--------------------------------------------------------------------------------AACGCTATGCCTGGGACACCATGCAGTACCAAATGGAACGTTAGTGAAAACAATCACGGATGACCAAATTGAAGTGACTAATGCTACTGAGCTGGTTCAAAGTTCCTCAACAGGTAGAATATGTAACAGTCCTCACCACATCCTTGATGGGAAAAATTGCACACTGATAGATGCTCTATTGGGAGACCCTCATTGTGATGACTTCCAAAACAAGGAATGGGACCTTTTTGTTGAACGAAGCACAGCCTACAGCAACTGCTACCCTTATTATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACCCTGGAATTCACCCAAGAAAGCTTCAATTGGACCGGAGTTACTCAAGATGGATCAAGCTATACTTGCAGAAGGAAATCTGTTAACAGTTTCTTTAGTAGATTAAATTGGTTGCATAATTTGGACTACAAATATCCAGCGCTGAACGTAACTATGCCAAACAATGACAAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACGGACAGGGACCAAACCAACCTATATGTTCAAGCATCAGGGAGAGTTACAGTCTCCACAAAAAGAAGCCAACAAACTGTAATCCCGAACATCGGATCTAGACCCTGGGTAAGGGGTGTCTCCAGCATAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGAAATCTAATTGCCCCTCGGGGTTACTTCAAAATACAAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAACTGCAATTCTGAATGCATTACTCCAAATGGAAGCATTCCCAATGACAAACCTTTTCAAAATGTAAACAGGATCACATATGGGGCCTGTCCAAGATATGTTAAGCAAAACACTCTGAAATTGGCAACAGGGATGCGGAATGTACCAGAGAAACAAACTAGAGGCATATTCGGCGCAATCGCAGGCTTCATAGAAAATGGTTGGGAGGGGATGGTGGACGGTTGGTACGGTTTCAGGCATCAAAATTCTGAAGGCACAGGACAAGCAGCAGATCTTAAAAGTACTCAAGCAGCAATCAACCAAATCACCGGGAAACTGAATAGAGTAATCAAGAAAACGAACGAGAAATTCCATCAAATCGAAAAAGAATTCTCAGAAGTAGAAGGGAGAATTCAGGACCTAGAGAAATACGTTGAAGACACTAAAATAGATCTCTGGTCTTACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATTTAACTGACTCAGAAATGAACAAACTGTTCGAAAGAACAAGGAAGCAACTGCGGGAAAATGCTGAGGACATGGGCAATGGTTGCTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATCAGAAATGGAACTTATGACCATGATGTATACAGAGACGAGGCATTAAACAATCGGTTCCAGATCAAAAGTGTTCAGCTGAAGTCAGGATACAAAGATTGGATCCTATGGATTTCCTTTGCCATGTCATGCTTTTTGCTTTGTGTTGTTCTGCTGGGGTTCATTATGTGGACCTGCCAAAAAGGCAACATTAAGTGCAACATTTGCATTTGA",
			"------------------------------------------------CAAAAACTTCCCGGAAATGACAACAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTACCAAACGGAACGATAGTGAAAACAATCACGAATGACCAAATTGAAGTTACTAATGCTACTGAGCTGGTTCAGAGTTCCTCAACAGGTGGAATATGCGACAGTCCTCATCAGATCCTTGATGGAGAAAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCAGTGTGATGGCTTCCAAAATAAGAAATGGGACCTTTTTGTTGAACGCAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACACTGGAGTTTAACAATGAAAGCTTCGATTGGACTGGAGTCACTCAGAATGGAACAAGCTCTGCTTGCAAAAGGAGATCTAATAAAAGTTTCTTTAGTAGATTGAATTGGTTGACCCACTTAAAATACAAATACCCAGCATTGAACGTGACTATGCCAAACAATGAAAAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACGGACAGTGACCAAATCAGCCTATATGCTCAAGCATCAGGAAGAATCACAGTCTCTACCAAAAGAAGCCAACAAACTGTAATCCCGAATATCGGATCTAGACCCAGGGTAAGGGATGTCTCCAGCCGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTCGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAAATGCAATTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCATTTCAAAATGTAAACAGGATCACATATGGGGCCTGTCCCAGATATGTTAAGCAAAACACTCTGAAATTGGCAACAGGGATGCGAAATGTACCAGAGAAACAAACTAGAGGCATATTTGGCGCAATCGCGGGTTTCATAGAAAATGGTTGGGAGGGAATGGTGGACGGTTGGTACGGTTTCAGGCATCAAAATTCTGAGGGCACAGGACAAGCAGCAGATCTCAAAAGCNCTCAAGCAGCAATGAAGACTATCATTG--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------",
			"ATGAAGACCATCATTGCTTTGAGCTACATTTTCTGTCTGGCTCTCGGCCAAGACCTTCCAGGAAATGACAACAGCACAGCAACGCTGTGCCTGGGACATCATGCGGTGCCAAACGGAACACTAGTGAAAACAATCACAGATGATCAGATTGAAGTGACTAATGCTACTGAGCTAGTTCAGAGCTCCTCAACGGGGAAAATATGCAACAATCCTCATCGAATCCTTGATGGAATAGACTGCACACTGATAGATGCTCTATTGGGGGACCCTCATTGTGATGTTTTTCAAAATGAGACATGGGACCTTTTCGTTGAACGCAGCAAAGCTTTCAGCAACTGTTACCCTTATGATGTGCCAGATTATGCCTCCCTTAGGTCACTAGTTGCCTCGTCAGGCACTCTGGAGTTTATCACTGAGGGTTTCACTTGGACTGGGGTCACTCAGAATGGGGGAAGCAATGCTTGCAAAAGGGGACCTGGTAGCGGTTTTTTCAGTAGACTGAACTGGTTGACCAAATCAGGAAGCACATATCCAGTGCTGAACGTGACTATGCCAAACAATGACAATTTTGACAAACTATACATTTGGGGGGTTCACCACCCGAGCACGAACCAAGAACAAACCAGCCTGTATGTTCAAGCATCAGGGAGAGTCACAGTCTCTACCAGGAGAAGCCGGCAAACTATAATCCCGAATATCGGGTCCAGACCCTGGGTAAGGGGTCTGTCTAGTAGAATAAGCATCTATTGGACAATAGTTAAGCCGGGAGACGTACTGGTAATTAATAGTAATGGGAACCTAATCGCTCCTCGGGGTTATTTCAAAATGCGCACTGGGAAAAGCTCAATAATGAGGTCAGATGCACCTATTGATACCTGTATTTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAGCCCTTTCAAAACGTAAACAAGATCACATATGGAGCATGCCCCAAGTATGTTAAGCAAAACACCCTGAAGTTGGCAACAGGGATGCGGAATGTACCAGAGAAACAAACTAGAGGCCTATTCGGCGCAATAGCAGGTTTCATAGAAAATGGTTGGGAGGGAATGATAGACGGTTGGTACGGTTTCAGGCATCAAAATTCTGAGGGCACAGGACAAGCAGCAGATCTTAAAAGCACTCAAGCAGCCATCGACCAAATCAATGGGAAATTGAACAGGGTAATCGAGAAGACGAACGAGAAATTCCATCAAATCGAAAAGGAATTCTCAGAAGTAGAAGGGAGAATTCAGGACCTCGAGAAATACGTTGAAGACACTAAAATAGATCTCTGGTCTTACAATGCGGAGCTTCTTGTCGCTCTGGAGAATCAACATACAATTGACCTGACTGACTCGGAAATGAACAAGCTGTTTGAAAAAACAAGGAGGCAACTGAGGGAAAATGCTGAAGACATGGGCAATGGTTGCTTCAAAATATACCACAAATGTGACAACGCTTGCATAGAGTCAATCAGAAATGGGACTTATGACCATAATGTATACAGAGACGAAGCATTAAACAACCGGTTTCAGATCAAAGGTGTTGAACTGAAGTCTGGATACAAAGACTGGATCCTGTGGATTTCCTTTGCCATATCATGCTTTTTGCTTTGTGTTGTTTTGCTGGGGTTCATCATGTGGGCCTGCCAGAGAGGCAACATTAGGTGCAACATTTGCATTTGA"
			]

		for reassortant_seq in reassortant_seqs:
			for v in self.viruses:
				dist = distance(Seq(reassortant_seq), v)
				if (dist < 0.02) and v.num_date>2005:
					remove_viruses.append(v)
					if self.verbose>1:
						print "\tremoving", v.strain

		self.viruses = MultipleSeqAlignment([v for v in self.viruses if v not in remove_viruses])'''

	def clean_outliers(self):
		"""Remove single outlying viruses"""
		remove_viruses = []
		outlier_strains = ["A/Sari/388/2006", "A/SaoPaulo/36178/2015", "A/Pennsylvania/40/2010", "A/Pennsylvania/14/2010", "A/Pennsylvania/09/2011", "A/OSAKA/31/2005", "A/Ohio/34/2012", "A/Kenya/170/2011", "A/Kenya/168/2011", "A/Indiana/21/2013", "A/Indiana/13/2012", "A/Indiana/11/2013", "A/Indiana/08/2012", "A/Indiana/06/2013", "A/India/6352/2012", "A/HuNan/01/2014", "A/Helsinki/942/2013", "A/Guam/AF2771/2011", "A/Chile/8266/2003", "A/Busan/15453/2009", "A/Nepal/142/2011", "A/Kenya/155/2011", "A/Guam/AF2771/2011"]

		for outlier_strain in outlier_strains:
			for v in self.viruses:
				if (v.strain == outlier_strain):
					remove_viruses.append(v)
					if self.verbose > 1:
						print "\tremoving", v.strain
		self.viruses = MultipleSeqAlignment([v for v in self.viruses if v not in remove_viruses])

	def clean(self):
		self.clean_generic()
		#self.clean_outbreaks()
		#print "Number of viruses after outbreak filtering:",len(self.viruses)
		#self.clean_reassortants()
		#print "Number of viruses after reassortant filtering:",len(self.viruses)
		#self.clean_outliers()
		#print "Number of viruses after outlier filtering:",len(self.viruses)

class H9_refine(tree_refine):
	def __init__(self, **kwargs):
		tree_refine.__init__(self, **kwargs)

		'''self.epitope_mask = ""
		if "epitope_masks_fname" in self.kwargs and "epitope_mask_version" in self.kwargs:
			epitope_map = {}
			with open(self.kwargs["epitope_masks_fname"]) as f:
				for line in f:
					(key, value) = line.split()
					epitope_map[key] = value
			if self.kwargs["epitope_mask_version"] in epitope_map:
				self.epitope_mask = epitope_map[self.kwargs["epitope_mask_version"]]
		self.epitope_mask = np.fromstring(self.epitope_mask, dtype='S1')				# epitope_mask is numpy array'''

	def refine(self):
		self.refine_generic()  # -> all nodes now have aa_seq, xvalue, yvalue, trunk, and basic virus properties
		self.add_H9_attributes()

	'''def epitope_sites(self, aa):
		aaa = np.fromstring(aa, 'S1')
		return ''.join(aaa[self.epitope_mask[:len(aa)]=='1'])

	def nonepitope_sites(self, aa):
		aaa = np.fromstring(aa, 'S1')
		return ''.join(aaa[self.epitope_mask[:len(aa)]=='0'])

	def receptor_binding_sites(self, aa):
		
		Receptor binding site mutations from Koel et al. 2014
		These are (145, 155, 156, 158, 159, 189, 193) in canonical HA numbering
		need to subtract one since python arrays start at 0
		
		sp = 16
		rbs = map(lambda x:x+sp-1, [145, 155, 156, 158, 159, 189, 193])					
		return ''.join([aa[pos] for pos in rbs])'''
	
	def get_total_peptide(self, node):	
		#the concatenation of signal peptide, HA1, HA1		
		return node.aa_seq['SigPep']+node.aa_seq['HA1']+node.aa_seq['HA2']

	'''def epitope_distance(self, aaA, aaB):
		Return distance of sequences aaA and aaB by comparing epitope sites
		epA = self.epitope_sites(aaA)
		epB = self.epitope_sites(aaB)
		distance = sum(a != b for a, b in izip(epA, epB))
		return distance

	def nonepitope_distance(self, aaA, aaB):
		Return distance of sequences aaA and aaB by comparing non-epitope sites
		neA = self.nonepitope_sites(aaA)
		neB = self.nonepitope_sites(aaB)
		distance = sum(a != b for a, b in izip(neA, neB))
		return distance'''

	#def receptor_binding_distance(self, aaA, aaB):
		#Return distance of sequences aaA and aaB by comparing receptor binding sites
		#neA = self.receptor_binding_sites(aaA)
		#neB = self.receptor_binding_sites(aaB)
		#distance = sum(a != b for a, b in izip(neA, neB))
		#return distance

	def add_H9_attributes(self):
		root = self.tree.seed_node
		#root_total_aa_seq = self.get_total_peptide(root)
		#for node in self.tree.postorder_node_iter():
			#total_aa_seq = self.get_total_peptide(node)
			#node.ep = self.epitope_distance(total_aa_seq, root_total_aa_seq)
			#node.ne = self.nonepitope_distance(total_aa_seq, root_total_aa_seq)
			#node.rb = self.receptor_binding_distance(total_aa_seq, root_total_aa_seq)


class H9_HI(HI_tree):
	def __init__(self, **kwargs):
		HI_tree.__init__(self, **kwargs)

class H9_fitness(fitness_model):
	def __init__(self, **kwargs):
		if 'predictors' in self.kwargs:
			predictor_input = self.kwargs['predictors']
			fitness_model.__init__(self, predictor_input = predictor_input, **kwargs)
		else:
			fitness_model.__init__(self, **kwargs)

	def annotate_fitness(self, estimate_frequencies = True):
		self.predict(estimate_frequencies=estimate_frequencies)


class H9_process(process, H9_filter, H9_clean, H9_refine, H9_fitness):
	"""docstring for H9_process, H9_filter"""
	def __init__(self,verbose = 0, force_include = None,
				force_include_all = False, max_global= True, **kwargs):
		self.force_include = force_include
		self.force_include_all = force_include_all
		self.max_global = max_global
		#self.max_globalh = max_globalh
		process.__init__(self, **kwargs)
		H9_filter.__init__(self,**kwargs)
		H9_clean.__init__(self,**kwargs)
		H9_refine.__init__(self,**kwargs)
		#H9_HI.__init__(self,**kwargs)
		H9_fitness.__init__(self,**kwargs)
		self.verbose = verbose

	def run(self, steps, viruses_per_month=50, raxml_time_limit = 1.0, lam_HI=2.0, lam_pot=0.3, lam_avi=2):
		if 'filter' in steps:
			print "--- Virus filtering at " + time.strftime("%H:%M:%S") + " ---"
			self.filter()
			'''if self.force_include is not None and os.path.isfile(self.force_include):
				with open(self.force_include) as infile:
					forced_strains = [fix_name(line.strip().split('\t')[0]).upper() for line in infile]
			else:
				forced_strains = []
			self.subsample(viruses_per_month,
				prioritize=forced_strains, all_priority=self.force_include_all,
				region_specific = self.max_global)
			self.add_older_vaccine_viruses(dt = 3)
			self.add_older_new_viruses(dt = 3)'''
			self.dump()
		else:
			self.load()
		if 'align' in steps:
			self.align()   	# -> self.viruses is an alignment object
		if 'clean' in steps:
			print "--- Clean at " + time.strftime("%H:%M:%S") + " ---"
			self.clean()   # -> every node as a numerical date
			self.dump()
		if 'tree' in steps:
			print "--- Tree	 infer at " + time.strftime("%H:%M:%S") + " ---"
			self.infer_tree(raxml_time_limit)  # -> self has a tree
			self.dump()
		if 'ancestral' in steps:
			print "--- Infer ancestral sequences " + time.strftime("%H:%M:%S") + " ---"
			self.infer_ancestral()  # -> every node has a sequence
			self.dump()
		if 'refine' in steps:
			print "--- Tree refine at " + time.strftime("%H:%M:%S") + " ---"
			self.refine()
			self.dump()
		if 'frequencies' in steps:
			print "--- Estimating frequencies at " + time.strftime("%H:%M:%S") + " ---"
			self.determine_variable_positions()
			#self.all_clade_frequencies(self, **kwargs)
			self.estimate_frequencies(tasks = ["mutations", "clades", "tree"])
			if 'genotype_frequencies' in steps:
					self.estimate_frequencies(tasks = ["genotypes"])
			self.dump()

		method = 'nnl1reg'
		if 'HI' in steps:
			print "--- Adding HI titers to the tree " + time.strftime("%H:%M:%S") + " ---"
			try:
				self.determine_variable_positions()
				self.map_HI(training_fraction=1.0, method = 'nnl1reg',
					lam_HI=lam_HI, lam_avi=lam_avi, lam_pot=lam_pot, map_to_tree=True)
				self.map_HI(training_fraction=1.0, method = 'nnl1reg', force_redo=True,
					lam_HI=lam_HI, lam_avi=lam_avi, lam_pot=lam_pot, map_to_tree=False)
			except:
				print("HI modeling failed!")
			#freqs = self.determine_HI_mutation_frequencies(threshold = 0.1)
			#self.frequencies["mutations"]["global"].update(freqs)
			self.dump()

		if 'fitness' in steps:
			print "--- Estimating fitnesses at " + time.strftime("%H:%M:%S") + " ---"
			self.annotate_fitness()
			self.dump()

		if 'export' in steps:
			#self.add_titers()
			self.temporal_regional_statistics()
			# exporting to json, including the H9 specific fields
			self.export_to_auspice(tree_fields = [
				'ep', 'ne', 'rb', 'aa_muts','accession','isolate_id', 'lab','db', 'country', 'dfreq', 'fitness', 'pred_distance',
				'dHI', 'cHI', 'mHI', 'mean_HI_titers', 'avidity_tree', 'avidity_mut', 'potency_mut', 'potency_tree', 'mean_potency_mut', 'mean_potency_tree', 'autologous_titers'],
				   annotations = ['Y439', 'G1', 'Ck-Bei', 'G9', 'Korea', 'Y280'])
			#self.export_fasta_alignment()
			#self.export_newick_tree()
			if params.html:
				self.generate_indexHTML()
			#self.export_HI_mutation_effects()
			#self.export_clade_frequencies()
			#self.export_viruses()

		if 'HIvalidate' in steps:
			from diagnostic_figures import tree_additivity_symmetry, fmts

			print "--- generating validation figures " + time.strftime("%H:%M:%S") + " ---"
			print "-- number of non-zero branch parameters: ",np.sum([n.dHI>1e-3 for n in self.tree.postorder_node_iter()])
			print "-- number of non-zero mutation parameters: ",np.sum([val>1e-3 for val in self.mutation_effects.values()])
			for model in ['tree', 'mutation']:
				try:
					tree_additivity_symmetry(self, model)
					for fmt in fmts: plt.savefig(self.htmlpath()+'HI_symmetry_'+model+fmt)
				except:
					print("Can't generate symmetry/additivity figures")
			try:
				self.slopes_muts = slope_vs_mutation(self)
			except:
				print("Couldn't derive slopes, probably to small time interval")
			self.generate_validation_figures(method)


if __name__=="__main__":
	all_steps = ['filter', 'align', 'clean', 'tree', 'ancestral', 'refine', 'frequencies', 'export']

	from process import parser
	import matplotlib.pyplot as plt
	plt.ion()
	params = parser.parse_args()

	lt = time.localtime()
	num_date = round(lt.tm_year+(lt.tm_yday-1.0)/365.0,2)
	params.time_interval = (num_date-params.years_back, num_date)
	if params.interval is not None and len(params.interval)==2 and params.interval[0]<params.interval[1]:
		params.time_interval = (params.interval[0], params.interval[1])
	dt= params.time_interval[1]-params.time_interval[0]
	params.pivots_per_year = 12.0 if dt<5 else 6.0
	steps = all_steps[all_steps.index(params.start):(all_steps.index(params.stop)+1)]
	if params.skip is not None:
		for tmp_step in params.skip:
			if tmp_step in steps:
				print "skipping",tmp_step
				steps.remove(tmp_step)

	# add all arguments to virus_config (possibly overriding)
	virus_config.update(params.__dict__)
	virus_config['serum_Kc'] = 0.003
	# pass all these arguments to the processor: will be passed down as kwargs through all classes
	myH9 = H9_process(**virus_config)
	if params.test:
		myH9.load()
	else:
		myH9.run(steps, viruses_per_month = virus_config['viruses_per_month'],
				   raxml_time_limit = virus_config['raxml_time_limit'],
				   #lam_HI = virus_config['lam_HI'],
				   #lam_avi = virus_config['lam_avi'],
				   #lam_pot = virus_config['lam_pot'],
				   )

