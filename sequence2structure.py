from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from Bio import PDB as pdb
import re
import csv
"""
sequence2structure.py

This program takes a given kinase domain sequence and determines its sequence
based on aligning multiple templates.

KNOWN ISSUES:
1) Paths are hard coded

TO DO:
1) Modify program so that information is automatically read in
2) Currently works in split up sections, but not all together -- need to fix this ASAP
"""

#Methods
"""
def get_codes(template):
	log.verbose()
	env = environ()
	env.io.atom_files_directory = './'
	#add chain information, just make it 'A'
	aln = alignment(env)
	for (code, chain) in ((template, 'A'), (template, 'A')):
	    mdl = model(env, file = code, model_segment = ('FIRST:' + chain, 'LAST:' + chain))
	    aln.append_model(mdl, atom_files = code, align_codes = code + chain)

def build_tree(protein):
	for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
	                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
	                                    ((1., 1., 1., 1., 1., 0.), True, False)):
	    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
	               rr_file='$(LIB)/as1.sim.mat', overhang=30,
	               gap_penalties_1d=(-450, -50),
	               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
	               dendrogram_file=protein_name+i+'.tree',
	               alignment_type='tree', # If 'progresive', the tree is not
	                                      # computed and all structues will be
	                                      # aligned sequentially to the first
	               feature_weights=weights, # For a multiple sequence alignment only
	                                        # the first feature needs to be non-zero
	               improve_alignment=True, fit=True, write_fit=write_fit,
	               write_whole_pdb=whole, output='ALIGNMENT QUALITY')
def salign(): 
	pdb_code = (pdb_name.split("-"[0])) 
	name = pdb_code[0] #changed this from hard coded 4F7S; it does not seem like this variable is used anywhere else
	chain = str(pdb_code[1])
	PIR = open('active.ali','w')
	PIR.write(">P1;{0}\n".format(pdb_name))
	PIR.write("structureX:{0}".format(header))
	PIR.write("{0}*\n\n".format(structure_sequence.strip()))
	PIR.write(">P1;{0}\n".format(protein_name))
	PIR.write("sequence:{0}".format(header))
	PIR.write("{0}*\n\n".format(full_sequence.strip()))
	PIR.close()

	aln.write(file=protein_name+i+'.pap', alignment_format='PAP') 
	aln.write(file=protein_name+i+'.ali', alignment_format='PIR')

	aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
	           rr_file='$(LIB)/as1.sim.mat', overhang=30,
	           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
	           gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
	           alignment_type='progressive', feature_weights=[0]*6,
	           improve_alignment=False, fit=False, write_fit=True,
	           write_whole_pdb=False, output='QUALITY')
"""
class PDB_info(object):
    
    #This class is used to assign meaning to specific elements in a given row of the .csv file
    
    def __init__(self, row):
        self.id = row[0] #id number of the pdb file
        self.protein = row[1] #protein name the pdb file is associated with
        self.complete = row[2] #yes or give missing residues
        self.conformation = row[3] #active or inactive?
        self.mutation = row[4] #is there a mutation?  If so, what are the details?

class Best_Template(object):
	def __init__(self,row):
		self.protein = row[0]
		self.template = row[1]

datafile = open('./JAK2_test.csv', 'r') #Opens the structures file for reading
datafile2 = open('./structures.csv', 'r')
datafile3 = open('./protein2template', 'r')

datareader = csv.reader(datafile2) #reads structures file
data = [] #initializes a list called data
for row in datareader:
    data.append(row) #adds an element to data for each row in structures.csv 

datareader3 = csv.reader(datafile3) #reads structures file
data3 = [] #initializes a list called data
for row in datareader3:
    data3.append(row) #adds an element to data for each row in structures.csv 

pdb_info = [PDB_info(item) for item in data]
template_info = [Best_Template(item) for item in data3]

"""
for i in range(len(template_info)):
	protein = template_info[i].protein
	template = template_info[i].template
	act = 0
	inact = 0
	for j in range(1, len(pdb_info)):
		if pdb_info[j].protein == template:
			if pdb_info[j].conformation == 'active':
				act = 1
			if pdb_info[j].conformation == 'inactive':
				inact = 1
	if act == 1 and inact == 1:
		print protein + ', good to go'
	elif act == 0 and inact == 1:
		print protein + ', missing active template'
	elif act == 1 and inact == 0:
		print protein + ', missing inactive template'
	elif act == 0 and inact == 0: 
		print protein + ', fuck'
"""

# need something that if "fuck," will search for the next closest sequence and use that as the template
# step 1: if "fuck," then look at the third highest 

clustal_lines = datafile.readlines()
index_line = clustal_lines[0]
index = index_line.split(',')
del index[0]
big_dict = {}
"""
def get_best_template():
	values_lines = clustal_lines[i]
	values = values_lines.split(',')
	protein_name = values[0]
	del values[0]
	num_values = [float(x) for x in values]
	#index = [protein_name + ':' + s for s in index]
	# put in hash to zip index and protein name
	small_dict = dict(zip(index, num_values))
	#big_dict.update((protein_name, small_dict))
	sorted_small_dict = sorted(small_dict, key=small_dict.get)
	act = 0
	inact = 0
	inact_struc = []
	act_struc = []
"""
for i in range(1, len(clustal_lines)):
	values_lines = clustal_lines[i]
	values = values_lines.split(',')
	protein_name = values[0]
	del values[0]
	num_values = [float(x) for x in values]
	#index = [protein_name + ':' + s for s in index]
	# put in hash to zip index and protein name
	small_dict = dict(zip(index, num_values))
	#big_dict.update((protein_name, small_dict))
	sorted_small_dict = sorted(small_dict, key=small_dict.get)
	print sorted_small_dict
	act = 0
	inact = 0
	inact_struc = []
	act_struc = []
	for j in range(1, 50):
		next_best_name = sorted_small_dict[-j]
		for k in range(0, 260):
			match = re.match(next_best_name, pdb_info[k].protein)
			if match != False and pdb_info[k].conformation == 'inactive':
				if inact != 1:
					inact_struc.append(next_best_name + '_inactive')
					inact = 1
					print 'yo'
				else: continue
			elif match != False and pdb_info[k].conformation == 'active':
				if act != 1:
					act_struc.append(next_best_name + '_active')
					act = 1
					print 'yo'
				else:
					continue
			else:
				continue
		if act == 1 and inact == 1:
			break
		else:
			continue

	print protein_name + ':' + next_best_name
	print small_dict[next_best_name]
	print inact_struc
	print act_struc

	# possible solution: a hash where the values are private hashes; each key in the hash is a protein; key = protein1 vs protein 2; value = sim score
	# editconf to fill in chain name
	# to verify --> make structure from sequence, compare to a full structure

	if act_struc[0] != protein_name:	
		#salign.py
		log.verbose()
		env = environ()
		env.io.atom_files_directory = './'
		#add chain information, just make it 'A'
		aln = alignment(env)
		for (code, chain) in ((next_best_name, 'A'), (next_best_name, 'A')):
		    mdl = model(env, file = code, model_segment = ('FIRST:' + chain, 'LAST:' + chain))
		    aln.append_model(mdl, atom_files = code, align_codes = code + chain)

		for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
		                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
		                                    ((1., 1., 1., 1., 1., 0.), True, False)):
		    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
		               rr_file='$(LIB)/as1.sim.mat', overhang=30,
		               gap_penalties_1d=(-450, -50),
		               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
		               dendrogram_file=protein_name+i+'.tree',
		               alignment_type='tree', # If 'progresive', the tree is not
		                                      # computed and all structues will be
		                                      # aligned sequentially to the first
		               feature_weights=weights, # For a multiple sequence alignment only
		                                        # the first feature needs to be non-zero
		               improve_alignment=True, fit=True, write_fit=write_fit,
		               write_whole_pdb=whole, output='ALIGNMENT QUALITY')

		pdb_code = (pdb_name.split("-"[0])) 
		name = pdb_code[0] #changed this from hard coded 4F7S; it does not seem like this variable is used anywhere else
		chain = str(pdb_code[1])
		PIR = open('active.ali','w')
		PIR.write(">P1;{0}\n".format(pdb_name))
		PIR.write("structureX:{0}".format(header))
		PIR.write("{0}*\n\n".format(structure_sequence.strip()))
		PIR.write(">P1;{0}\n".format(protein_name))
		PIR.write("sequence:{0}".format(header))
		PIR.write("{0}*\n\n".format(full_sequence.strip()))
		PIR.close()

		aln.write(file=protein_name+i+'.pap', alignment_format='PAP') 
		aln.write(file=protein_name+i+'.ali', alignment_format='PIR')

		aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
		           rr_file='$(LIB)/as1.sim.mat', overhang=30,
		           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
		           gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
		           alignment_type='progressive', feature_weights=[0]*6,
		           improve_alignment=False, fit=False, write_fit=True,
		           write_whole_pdb=False, output='QUALITY')

		#align2d_mult.py
		log.verbose()
		env = environ()

		env.libs.topology.read(file='$(LIB)/top_heav.lib')

		# Read aligned structure(s):
		aln = alignment(env)
		aln.append(file=protein_name+i+'.ali', align_codes='all')
		aln_block = len(aln)

		# Read aligned sequence(s):
		aln.append(file=protein_name + '.ali', align_codes=protein_name)

		# Structure sensitive variable gap penalty sequence-sequence alignment:
		aln.salign(output='', max_gap_length=20,
		           gap_function=True,   # to use structure-dependent gap penalty
		           alignment_type='PAIRWISE', align_block=aln_block,
		           feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
		           gap_penalties_1d=(-450, 0),
		           gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
		           similarity_flag=True)

		aln.write(file=protein_name+'-mult.ali', alignment_format='PIR')
		aln.write(file=protein_name+'-mult.pap', alignment_format='PAP')

		#model_mult.py
		env = environ()
		a = automodel(env, alnfile=protein_name+'-mult.ali',
		              knowns=(next_best_name), sequence=protein_name)
		a.starting_model = 1
		a.ending_model = 5
		a.make()

		#evaluate_model.py
		log.verbose()    # request verbose output
		env = environ()
		env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
		env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

		# read model file
		mdl = complete_pdb(env, protein_name+'.B99990001.pdb')

		# Assess all atoms with DOPE:
		s = selection(mdl)
		s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=protein_name+'.profile',
		              normalize_profile=True, smoothing_window=15)
	elif inact_struc != protein_name:
		continue
	else: continue
