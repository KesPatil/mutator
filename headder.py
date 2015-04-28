#!/usr/bin/env 
from modeller import *
from modeller.optimizers import molecular_dynamics, conjugate_gradients
from modeller.automodel import *
from modeller.scripts import complete_pdb
from Bio import PDB as pdb
import re
import csv
import os, glob
import shutil

"""
headder.py

header adder

Adds a new header to completed structure files by putting in the missing residue and its corresponding number in the sequence.
"""

#function to add lines to top of the file
def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

#class to specify each element in CSV
class PDB_info(object):
    """
    This class is used to assign meaning to specific elements in a given row of the .csv file
    """
    def __init__(self, row):
        self.id = row[0] #id number of the pdb file
        self.protein = row[1] #protein name the pdb file is associated with
        self.complete = row[2] #yes or give missing residues
        self.conformation = row[3] #active or inactive?
        self.mutation = row[4] #is there a mutation?  If so, what are the details?

#reads in info from csv file
datafile = open('./temp.csv', 'r') #Opens the structures file for reading
datareader = csv.reader(datafile) #reads structures file
data = [] #initializes a list called data
for row in datareader:
    data.append(row) #adds an element to data for each row in structures.csv 

#parses csv data using PDB_info class
pdb_info = [PDB_info(item) for item in data]

for i in range (1, len(pdb_info)):
	#assigns variable names to pdb_info elements
	pdb_name = pdb_info[i].id #saves given pdb name as a variable
	protein_name = pdb_info[i].protein #saves given protein name as a variable
	complete = pdb_info[i].complete #saves yes or no for complete
	structure_conf = pdb_info[i].conformation #saves active or inactive for conformation
	mutation = pdb_info[i].mutation

    #gives location of the pdb file
	pdb_file = './PDBs/'+pdb_name+'.pdb'


    #gets name of the structure file
	if structure_conf == 'active':
		if complete == 'yes':
			if mutation == 'no_mutations':
				structure_file = './actives/complete/'+protein_name+'_active.pdb'
			else:
				different_mutations = re.split('&', mutation)
				structure_file = './actives/complete/'+protein_name+'_active'+'.pdb'
		else:
			if mutation == 'no_mutations':
				structure_file = './actives/incomplete/'+protein_name+'_active.pdb'
			else:
				structure_file = './actives/incomplete/'+protein_name+'_active'+'.pdb'
	else:
		if complete == 'yes':
			if mutation == 'no_mutations':
				structure_file = './inactives/complete/'+protein_name+'_inactive.pdb'
			else:
				structure_file = './inactives/complete/'+protein_name+'_inactive'+'.pdb'
		else:
			if mutation == 'no_mutations':
				structure_file = './inactives/incomplete/'+protein_name+'_inactive.pdb'
			else:
				structure_file = './inactives/incomplete/'+protein_name+'_inactive'+'.pdb'


    #regular expression for the header (COMPLETE THIS)
	fp = open(pdb_file)
	lines = fp.readlines()
	good_lines = []
	for n in range(0, 10):
		if re.match('HEADER', lines[n]) is not None:
			#get rid of new line character
			better_line = re.sub('\\n', '', lines[n])
			good_lines.append(better_line)
		elif re.match('REMARK 300', lines[n]) is not None:
			#get rid of new line character
			better_line = re.sub('\\n', '', lines[n])
			good_lines.append(better_line)
		else: break
	

	rev_good_lines = []
	for i in reversed(good_lines):
		rev_good_lines.append(i)

	print rev_good_lines

    #concat the header and structure file (COMPLETE THIS)
	for m in range(0, len(rev_good_lines)):
		line_prepender(structure_file, rev_good_lines[m])

    #get missing residue ranges


    #parse sequence into counter (res #, res)

    #print missing res w/ corresponding number and concat

