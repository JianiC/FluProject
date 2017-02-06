import matplotlib as mpl
mpl.use('pdf')
import time, re, os
from virus_filter import flu_filter, fix_name
from virus_clean import virus_clean
from tree_refine import tree_refine
from tree_titer import HI_tree
from fitness_model import fitness_model
from H9_process import H9_refine as H4_refine
from process import process, virus_config
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import numpy as np
from itertools import izip

# numbering starting at methionine including the signal peptide
'''
sp = 17
epitope_mask = np.array(['1' if pos in [141,142,145,146,172,176,178,179,180,181,183,184,185, #Sa
										170,173,174,177,206,207,210,211,212,214,216,		 #Sb
										183,187,191,196,221,225,254,258,288,				 #Ca1
										154,157,158,159,161,163,238,239,242,243,			 #Ca2
										87, 88, 90, 91, 92, 95, 96, 98, 99, 100, 132, 139	 #Cb
									   ]
						else '0' for pos in xrange(1,1725)])

receptor_binding_sites = [x-1 for x in [159,169,170,172,173,203,207]]
'''

virus_config.update({
	# data source and sequence parsing/cleaning/processing
	'virus':'H4',
	'alignment_file':'/Users/yujiazhou/Documents/nextflu/H9_nextflu-master/augur/src/data/H4_gisaid_epiflu_sequence.fasta',
	'outgroup':'A/Duck/Czechoslovakia/1956',
	#'force_include':'/Users/yujiazhou/Documents/nextflu/H9_nextflu-master/augur/src/data/H4_HI_strains.txt',
	'force_include_all':True,
	'date_spec':'year',
	'max_global':True,   # sample as evenly as possible from different geographic regions

	'cds':[0,None], # define the HA start i n 0 numbering

	# define relevant clades in canonical HA1 numbering (+1)
	"""
	'clade_designations': {
		'2': [('HA1', 125, 'N'), ('HA1', 134 ,'A'), ('HA1', 183, 'S'), ('HA1', 31,'D'), ('HA1', 172,'N'), ('HA1', 186,'T')],
		'3': [('HA1', 134 ,'T'), ('HA1', 183, 'P')],
		'4': [('HA1', 125, 'D'), ('HA1', 134 ,'A'), ('HA1', 183, 'S')],
		'5': [('HA1', 87, 'N'), ('HA1', 205, 'K'), ('HA1', 216, 'V'), ('HA1', 149, 'L')],
		'6': [('HA1', 185,'T'),  ('HA1', 97, 'N'), ('HA1', 197, 'A')],
		'6c':[('HA1', 234,'I'),  ('HA1', 97, 'N'), ('HA1', 197, 'A'), ('HA1', 283,'E')],
		'6b':[('HA1', 163,'Q'),  ('HA1', 256, 'T'), ('HA1', 197, 'A'), ('HA1', 283,'E')],
		'7': [('HA1', 143,'G'),  ('HA1', 97, 'D'), ('HA1', 197, 'T')],
		'8': [('HA1', 186,'T'),  ('HA1', 272,'A')],
		'84N':[('HA1', 163,'Q'),  ('HA1', 256, 'T'), ('HA1', 197, 'A'), ('HA1', 283,'E'), ('SigPep', 13, 'T'), ('HA1', 84, 'N')]
		},
	"""
	#'HI_fname':'/Users/yujiazhou/Documents/nextflu/H9_nextflu-master/augur/src/data/H4_HI_titers.txt',
	'auspice_prefix':'H4_',
	'html_vars': {'coloring': 'region, date',
				  'gtplaceholder': 'HA1 positions...',
				  'freqdefault': ''},
	'js_vars': {'LBItau': 0.0005, 'LBItime_window': 0.5, 'dfreq_dn':2},
	'layout':'auspice',
	})


class H4_filter(flu_filter):
	def __init__(self,min_length = 0, **kwargs):
		'''
		parameters
		min_length  -- minimal length for a sequence to be acceptable
		'''
		flu_filter.__init__(self, **kwargs)
		self.min_length = min_length
		self.vaccine_strains =[
			{
				'strain':		'A/ruddyshelduck/Mongolia/1626/2010',
				'isolate_id':	'EPI-ISL-149888',
				'date':			'2010-09-11',
				'lab':			'Other Database Import',
				'region':		'NorthAsia',
				'seq':			'ATGCTATCAATTGTGATTTTGTTTCTGCTTGTTGCAGAGAGCTCTTCTCAAAACTACACAGGAAACCCTGTGATATGCATGGGACATCATGCTGTGGCCAATGGGACTATGGTAAAGACCCTTACTGATGATCAAGTGGAAGTGGTCACTGCACAAGAATTGGTGGAATCACAGAACCTCCCGGAACTATGCCCGAGTCCTCTAAGACTAGTCGATGGCCAGACCTGTGATATCATCAATGGAGCCTTAGGAAGCCCAGGATGTGACCATTTGAATGGTGCTGAATGGGACATTTTCATAGAAAGGCCCAATGCAGTGGACACTTGCTATCCATTTGATGTGCCAGATTATCAGAGCCTAAGGAGCATACTCGCCAACAATGGGAAATTCGAATTCATTGCTGAAGAATTCCAATGGAGCACCGTGAAGCAAAATGGCAAGTCCGGGGCCTGCAAGAGGGCAAATGTGAACGATTTCTTTAATAAACTGAATTGGCTCGTGAAGTCAGACGGGAATGCATACCCTCTCCAGAATTTGACAAAAGTAAACAACGGTGATTACGCGAGGCTTTACATCTGGGGAGTTCACCACCCTTCGACGGATACCGAGCAAACCGATCTGTACAAGAACAATCCTGGTAGGGTCACTGTATCTACCAAAATCAGTCAAACAAGTGTAGTGCCCAACATTGGCAGCAGACCTTGGGTGAGAGGACAAAGTGGCAGAATCAGCTTCTATTGGACTATTGTAGAGCCTGGAGATTTGATAGTCTTCAACACAATAGGAAATTTAATTGCCCCAAGAGGACATTACAAATTAAACAGTCAGAAGAAGAGCACAATTCTGAACACTGCGACTCCCATAGGCTCATGTGTCAGTAAATGTCATACAGACAAAGGTTCTCTCTCTACCACCAAGCCCTTTCAAAATATCTCAAGGATAGCAGTTGGAGATTGTCCCAAATATGTTAAACAAGGCTCCCTAAAACTTGCAACTGGGATGAGAAATATCCCTGAAAAGGCATCAAGAGGGCTTTTTGGGGCAATAGCTGGGTTCATAGAGAATGGATGGCAAGGTCTGATTGATGGTTGGTATGGCTTCAGACACCAAAATGCAGAAGGAACAGGAACAGCTGTTGATCTAAAATCCACTCAGGCAGCCATCGATCAAATCAATGGAAAACTCAATCGTCTTATTGAGAAAACAAACGAGAAATACCATCAAATCGAAAAAGAATTCGAACAAGTTGAAGGAAGAATCCAAGACCTGGAGAAGTATGTTGAAGACACAAAGATTGATCTATGGTCATATAATGCAGAGCTATTAGTCGCTCTGGAAAACCAGCATACTATAGATGTGACTGACTCGGAGATGAACAAGCTCTTTGAAAGAGTAAGGCGACAACTCAGGGAGAATGCTGAAGACAGAGGAAATGGGTGTTTTGAAATATTCCACAAATGTGACAACAACTGCATTGAAAGCATTCGGAATGGGACCTATGATCATGATGTTTATAGAGATGAAGCGATCAACAATCGATTCCAAATACAGGGAGTCAAATTGACCCAGGGATACAAGGACATCATCCTTTGGATTTCGTTCTCCATATCATGCTTTTTGCTCGTAGCACTGCTTTTGGCCTTCATTTTGTGGGCTTGTCAGAACGGAAACATCCGGTGCCAGATT---TGCATTTGA',
			}
		]
		tmp_outgroup = SeqIO.read('/Users/yujiazhou/Documents/nextflu/H9_nextflu-master/augur/source-data/H4_outgroup.gb', 'genbank')
		genome_annotation = tmp_outgroup.features
		self.cds = {x.qualifiers['gene'][0]:x for x in genome_annotation
				if 'gene' in x.qualifiers and x.type=='CDS' and
				x.qualifiers['gene'][0] in ['SigPep', 'HA1', 'HA2']}
		self.outgroup = {
			'strain': 'A/Duck/Czechoslovakia/1956',
			'db': 'OtherDatabaseImport',
			'accession': 'EPI-ISL-70104',
			'date': '1956-01-01',
			'country': 'Czech Republic',
			'region': 'Europe',
			'seq': str(tmp_outgroup.seq).upper()
		}


class H4_clean(virus_clean):
	def __init__(self,**kwargs):
		virus_clean.__init__(self, **kwargs)

	def clean_outbreaks(self):
		"""Remove duplicate strains, where the geographic location, date of sampling and sequence are identical"""
		virus_hashes = set()
		new_viruses = []
		for v in self.viruses:
			geo = re.search(r'A/([^/]+)/([^/]+)/', v.strain).group(2)
			if geo:
				vhash = (geo, v.date, str(v.seq))
				if vhash not in virus_hashes:
					new_viruses.append(v)
					virus_hashes.add(vhash)

		self.viruses = MultipleSeqAlignment(new_viruses)
		return new_viruses

	def clean(self):
		self.clean_generic()
		#self.clean_outbreaks()
		#print "Number of viruses after outbreak filtering:",len(self.viruses)
		#self.clean_outliers()
		#self.clean_outlier_strains()
		#print "Number of viruses after outlier filtering:",len(self.viruses)

class H4_process(process, H4_filter, H4_clean, H4_refine, fitness_model):
	"""docstring for H4_process, H4_filter"""
	def __init__(self,verbose = 0, force_include = None,
				force_include_all = False, max_global= True, **kwargs):
		self.force_include = force_include
		self.force_include_all = force_include_all
		self.max_global = max_global
		process.__init__(self, **kwargs)
		H4_filter.__init__(self,**kwargs)
		H4_clean.__init__(self,**kwargs)
		H4_refine.__init__(self,**kwargs)
		#HI_tree.__init__(self,**kwargs)
		fitness_model.__init__(self,**kwargs)
		self.verbose = verbose

	def run(self, steps, viruses_per_month=50, raxml_time_limit=1.0, lam_HI=2.0, lam_pot=0.3, lam_avi=2.0):
		if 'filter' in steps:
			print "--- Virus filtering at " + time.strftime("%H:%M:%S") + " ---"
			self.filter()
			if self.force_include is not None and os.path.isfile(self.force_include):
				with open(self.force_include) as infile:
					forced_strains = [fix_name(line.strip().split('\t')[0]).upper() for line in infile]
			else:
				forced_strains = []
			self.subsample(viruses_per_month,
				prioritize=forced_strains, all_priority=self.force_include_all,
				region_specific = self.max_global)
			self.add_older_vaccine_viruses(dt = 6)
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
			self.estimate_frequencies(tasks = ["mutations", "clades", "tree"])
			if 'genotype_frequencies' in steps:
					self.estimate_frequencies(tasks = ["genotypes"])
			self.dump()
		if 'export' in steps:
			#self.add_titers()
			self.temporal_regional_statistics()
			# exporting to json, including the H4 specific fields
			self.export_to_auspice(tree_fields = [
				'aa_muts','accession','isolate_id', 'lab','db', 'country',
				'avidity_tree','avidity_mut', 'potency_mut', 'potency_tree', 'mean_potency_mut', 'mean_potency_tree', 'autologous_titers'],
                   annotations = ['5', '6', '6b', '6c', '7', '84N'])
			if params.html:
				self.generate_indexHTML()
			#self.export_HI_mutation_effects()

if __name__=="__main__":
	all_steps = ['filter', 'align', 'clean', 'tree', 'ancestral', 'refine', 'frequencies', 'export']
	from process import parser
	params = parser.parse_args()

	lt = time.localtime()
	num_date = round(lt.tm_year+(lt.tm_yday-1.0)/365.0,2)
	print "num_date", num_date
	params.time_interval = (num_date-params.years_back, num_date)
	print "params.time_interval", params.time_interval
	if params.interval is not None and len(params.interval)==2 and params.interval[0]<params.interval[1]:
		params.time_interval = (params.interval[0], params.interval[1])
	dt= params.time_interval[1]-params.time_interval[0]
	params.pivots_per_year = 12.0 if dt<5 else 6.0 if dt<10 else 3.0
	steps = all_steps[all_steps.index(params.start):(all_steps.index(params.stop)+1)]
	if params.skip is not None:
		for tmp_step in params.skip:
			if tmp_step in steps:
				print "skipping",tmp_step
				steps.remove(tmp_step)


	# add all arguments to virus_config (possibly overriding)
	virus_config.update(params.__dict__)
	# pass all these arguments to the processor: will be passed down as kwargs through all classes
	myH4 = H4_process(**virus_config)
	if params.test:
		myH4.load()
	else:
		myH4.run(steps, viruses_per_month = virus_config['viruses_per_month'],
				raxml_time_limit = virus_config['raxml_time_limit'],
				   #lam_HI = virus_config['lam_HI'],
				   #lam_avi = virus_config['lam_avi'],
				   #lam_pot = virus_config['lam_pot'],
				   )

