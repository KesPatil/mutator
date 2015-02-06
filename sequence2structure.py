from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
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

datafile = open('./kinase_clustalo.csv', 'r') #Opens the structures file for reading
"""
datareader = csv.reader(datafile) #reads structures file
data = [] #initializes a list called data
for row in datareader:
    data.append(row) #adds an element to data for each row in structures.csv 
"""

def second_largest(numbers):
    count = 0
    m1 = m2 = float('-inf')
    for x in numbers:
        count += 1
        if x > m2:
            if x >= m1:
                m1, m2 = x, m1            
            else:
                m2 = x
    return m2 if count >= 2 else None

clustal_lines = datafile.readlines()
index_line = clustal_lines[0]
index = index_line.split(',')

for i in range(1, len(clustal_lines)):
	values_lines = clustal_lines[i]
	values = values_lines.split(',')
	protein_name = values[0]
	
	num_values = []
	for j in range (1, len(values)):
		num_values.append(float(values[j]))
	#for k in range(0, len(int_values)):
		#if int_values[k] == 100:
		#	int_values[k] = 0
		#if int_values[k] == max(int_values):
		#	best_match_index = k
	best_match = second_largest(num_values)
	best_match_index = num_values.index(best_match) + 1
	best_template = index[best_match_index]
	
	#best_match = int_values[best_match_index]
	print protein_name, best_match, best_match_index, best_template

	
#	best_match = second_largest(int_values)
#	re.search(best_match, int_values)
#	print int_values
	

"""
#salign.py
log.verbose()
env = environ()
env.io.atom_files_directory = './'

aln = alignment(env)
for (code, chain) in (('2mdh', 'A'), ('1bdm', 'A'), ('1b8p', 'A')):
    mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(mdl, atom_files=code, align_codes=code+chain)

for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50),
               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
               dendrogram_file='fm00495.tree',
               alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
               feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
               improve_alignment=True, fit=True, write_fit=write_fit,
               write_whole_pdb=whole, output='ALIGNMENT QUALITY')

aln.write(file='fm00495.pap', alignment_format='PAP')
aln.write(file='fm00495.ali', alignment_format='PIR')

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
aln.append(file='fm00495.ali', align_codes='all')
aln_block = len(aln)

# Read aligned sequence(s):
aln.append(file='TvLDH.ali', align_codes='TvLDH')

# Structure sensitive variable gap penalty sequence-sequence alignment:
aln.salign(output='', max_gap_length=20,
           gap_function=True,   # to use structure-dependent gap penalty
           alignment_type='PAIRWISE', align_block=aln_block,
           feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
           gap_penalties_1d=(-450, 0),
           gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
           similarity_flag=True)

aln.write(file='TvLDH-mult.ali', alignment_format='PIR')
aln.write(file='TvLDH-mult.pap', alignment_format='PAP')

#model_mult.py
env = environ()
a = automodel(env, alnfile='TvLDH-mult.ali',
              knowns=('1bdmA','2mdhA','1b8pA'), sequence='TvLDH')
a.starting_model = 1
a.ending_model = 5
a.make()

#evaluate_model.py
log.verbose()    # request verbose output
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# read model file
mdl = complete_pdb(env, 'TvLDH.B99990001.pdb')

# Assess all atoms with DOPE:
s = selection(mdl)
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='TvLDH.profile',
              normalize_profile=True, smoothing_window=15)
"""
