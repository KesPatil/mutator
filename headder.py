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
	ppb = pdb.PPBuilder() #peptide class to get sequence
	last = 10000

    #gives location of the pdb file
	pdb_file = './PDBs/'+pdb_name+'.pdb'
	parser = pdb.PDBParser()
	struct = parser.get_structure("name",pdb_file) #read in pdb file using PDBParser

    #gets name of the structure file
	if structure_conf == 'active':
		if complete == 'yes':
			structure_file = './actives/complete/'+protein_name+'_active.pdb'
		else:
			structure_file = './actives/incomplete/'+protein_name+'_active.pdb'
	else:
		if complete == 'yes':
			structure_file = './inactives/complete/'+protein_name+'_inactive.pdb'
		else:
			structure_file = './inactives/incomplete/'+protein_name+'_inactive.pdb'


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

    #get missing residue ranges

	structure_sequence = ''
	first_range = []
	last_range = []
	for seq in ppb.build_peptides(struct):
	    print seq.get_sequence()
	#read in the full sequence from the pdb file
	full_sequence = ''
	first = lines[0]
	header = re.split('HEADER\W+',first)[1] #uses a modified version of PDB file
	header_list = header.split(':')
	first_res = int(header_list[1])
	last_res = int(header_list[3])

	for seq in ppb.build_peptides(struct):
	    #use this re to get the chain breaks
	    search = re.search('start=([0-9]{1,5}).+end=([0-9]{1,5})',"{0}".format(seq))
	    first_range.append(search.groups()[0])
	    last_range.append(search.groups()[1])
	    first = search.groups()[0]
	    diff = int(first)-int(last)-1
	    if(diff > 0): #put in dashes for missing residues
	        structure_sequence += diff*'-'
	    last = search.groups()[1]
	    structure_sequence += seq.get_sequence()

	    #print (int(first)-21),(int(last)-21)
	print structure_sequence

	first_range = map(int, first_range) #makes this an integer array
	last_range = map(int, last_range) #makes this an integer array
	first_res_in_range = first_range.pop(0) #gets rid of the first element
	last_res_in_range = last_range.pop(-1) #gets rid of the last element
	first_missing = [x + 1 for x in last_range] #will use this to make missing residue ranges
	last_missing = [x - 1 for x in first_range] #will use this to make missing residue ranges

	for i in range(0, len(first_missing)):
		print first_missing[i], last_missing[i] + 1

    #parse sequence into counter (res #, res)
	for index in range(1,10):
		split_line = re.split('REMARK 300 ',lines[index])
		if split_line[0] == '':
			full_sequence += split_line[1]
	full_sequence.rstrip('\n')
	print full_sequence
	full_sequence_list = list(full_sequence)
	mis_res_list = []
	for i in range(0, len(first_missing)):
		start = first_missing[i] - first_res
		end = last_missing[i] + 1 - first_res
		numbers = range(start, end)
		print numbers
		for j in range(start, end):
			real_start = j - start
			mis_res_list.append(str(numbers[real_start])+': '+full_sequence_list[j])
		print mis_res_list
	joiner = ', '.join(mis_res_list)
	line_prepender(structure_file, mutation)
	line_prepender(structure_file, 'Mutations:')
	line_prepender(structure_file, joiner)
	line_prepender(structure_file, 'Missing Residues:')
	for m in range(0, len(rev_good_lines)):
		line_prepender(structure_file, rev_good_lines[m])	
    #print missing res w/ corresponding number and concat

