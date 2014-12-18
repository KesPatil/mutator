#!/usr/bin/env python
from Bio import PDB as pdb
import os
import random
import csv
import re
"""
1) Purpose of this program: delete 10 residues from random places in complete structures
2) Build back the structures
3) Write a .pml file in order to compare RMSDs
"""
datafile = open('./active_verifier.csv', 'r') #Opens the structures file for reading
datareader = csv.reader(datafile) #reads structures file
data = [] #initializes a list called data
for row in datareader:
    data.append(row) #adds an element to data for each row in structures.csv 

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
pdb_info = [PDB_info(item) for item in data]
for i in range(1, len(pdb_info)):
	pdb_name = pdb_info[i].id #saves given pdb name as a variable
	protein_name = pdb_info[i].protein #saves given protein name as a variable
	complete = pdb_info[i].complete #saves yes or no for complete
	structure_conf = pdb_info[i].conformation #saves active or inactive for conformation
	mutation = pdb_info[i].mutation
	print pdb_name
	pdb_file = './PDBs/'+pdb_name+'.pdb'
	fp = open(pdb_file)
	lines = fp.readlines()
	first = lines[0]
	header = re.split('HEADER\W+',first)[1] #uses a modified version of PDB file
	header_list = header.split(':')
	first_res = int(header_list[1])
	last_res = int(header_list[3])
	fp.close()
	fp = open(pdb_file)
	regex = re.split('\n'+'ATOM', fp.read())
	adder = regex[0]
	adder = re.sub(pdb_name, 'cropped_'+protein_name+'_active', adder)
	random_int = random.randint(first_res, last_res - 21)
	residues_to_remove = []
	for j in range(1, 5):
		residues_to_remove.append(random_int + j)
	parser = pdb.PDBParser()
	structure = parser.get_structure("name", pdb_file)
	first_model = structure[0]
	for chain in first_model:
	    for id in residues_to_remove:
	        chain.detach_child((' ', id, ' '))
	io = pdb.PDBIO()
	io.set_structure(first_model)
	io.save('./cropped/cropped_'+protein_name+'_active.pdb')
	template = open('./cropped/cropped_'+protein_name+'_active.pdb', 'r')
	saver = template.read()
	template.close()
	template = open('./cropped/cropped_'+protein_name+'_active.pdb', 'w')
	template.write(adder + '\n')
	template.write(saver)
	template_briefname = ('cropped_'+protein_name)
	template_csv = open('cropped_actives.csv','a')
	template_csv.write(template_briefname+','+template_briefname+'_modelled,'+'no,'+structure_conf+','+mutation+"\n")
	
	"""
	class Prepender:

	    def __init__(self, fname, mode='w'):
	        self.__write_queue = []
	        self.__f = open(fname, mode)

	    def write(self, s):
	        self.__write_queue.insert(0, s)

	    def close(self):
	        self.__exit__(None, None, None)

	    def __enter__(self):
	        return self

	    def __exit__(self, type, value, traceback):
	        if self.__write_queue: 
	            self.__f.writelines(self.__write_queue)
	        self.__f.close()

	with Prepender(template) as f:
	    f.write(first_model)
	    f.write(adder+'\n')
	"""