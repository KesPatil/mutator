#!/usr/bin/env 
from modeller import *
from modeller.optimizers import molecular_dynamics, conjugate_gradients
from modeller.automodel import *
from modeller.scripts import complete_pdb
from Bio import PDB as pdb
import re
import csv
import os, glob
import gromacs
import random
import numpy

"""
loop_builder.py

This program fills in missing residues from (correctly formatted) pdb files.
It reads in the coordinates and detects any missing residues and then
gives an alignment file that can be used for loop modeling. It currently
does not do a very good job, as modeller is both poorly executed and also
I do not have the practical knowlede to make it do exactly what I want.

Known issues: 
    the paths are hard coded, you will need to modify them
    tested 3CQU-A and 3O96-A, got RMSD of 1.5 -- not as good as could be
    when the first residue is missing, it can't write the alignment file

To do:
    fix issues with knots!
        *** (DONE) ***
        *** Solution: make sure the dashes in the "active.ali" file are of
        correct length ***
    modify way it reads in sequence so that models are automatically given
        the correct name including active/inactive designation *** (DONE) ***
    put in loop to automatically put in regions to modify (now hard coded)
        *** (DONE) ***
    put it in a loop so that it reads the active/inactive templates and then
        create models in a separate directory *** (DONE) ***
        *** (DONE) ***
    put in a sanity check to make sure there are no mutations in the pdb 
        structures (we know there are some) and use modeller to mutate
        these back 
        *** how to call mutate_model.py? ***
        *** should I make a new text file that will save the mutated models,
        then feed that text file into mutate_model.py? ***
    maybe put in a sanity check at the end that compares models of similar
        sequence or structure
    maybe autoload into pymol?
    maybe a sanity check for whether protein_name is correct?
    sanity check to make sure that it is actually complete or incomplete
"""

datafile = open('./temp.csv', 'r') #Opens the structures file for reading
datareader = csv.reader(datafile) #reads structures file
data = [] #initializes a list called data
for row in datareader:
    data.append(row) #adds an element to data for each row in structures.csv 

datafile2 = open('./kinase_domain_info_uniprot.csv', 'r') #Opens the structures file for reading
datareader2 = csv.reader(datafile2) #reads structures file
data2 = [] #initializes a list called data
for row in datareader2:
    data2.append(row) #adds an element to data for each row in structures.csv 

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

class Kinase_domain_info(object):
    """
    This class is used to assign meaning to specific elements in a given row of the .csv file
    """
    def __init__(self, row):
        self.protein_name = row[0] 
        self.kdstart = row[1] 
        self.kdend = row[2] 
        self.ploopstart = row[3] 
        self.ploopend = row[4] 
        self.alphacstart = row[5]
        self.alphacend = row[6]
        self.activationloopstart = row[7]
        self.activationloopend = row[8]

pdb_info = [PDB_info(item) for item in data]
kinase_domain_info = [Kinase_domain_info(item) for item in data2] 
#print pdb_info[0].id

for i in range(1, len(pdb_info)):
    pdb_name = pdb_info[i].id #saves given pdb name as a variable
    protein_name = pdb_info[i].protein #saves given protein name as a variable
    complete = pdb_info[i].complete #saves yes or no for complete
    structure_conf = pdb_info[i].conformation #saves active or inactive for conformation
    mutation = pdb_info[i].mutation

    """
    search for a given string where protein_name == kinase_domain_info.protein_name
    then go to that line 
    """

    print pdb_name
    pdb_file = ('./PDBs/'+pdb_name+'.pdb')
    if os.path.isfile(pdb_file) != True: #if there is no pdb_file, then make a note of it and continue with next pdb
        missing = open('missing_PDBs.csv','a')
        missing.write(pdb_name+','+protein_name+','+complete+','+structure_conf+','+mutation+"\n")
    else:
        fp = open(pdb_file)
        parser = pdb.PDBParser()
        struct = parser.get_structure("name",pdb_file) #read in pdb file using PDBParser
        ppb = pdb.PPBuilder() #peptide class to get sequence
        last = 100000 #make sure first iter through loop has no dashes -
        structure_sequence = ''
        first_range = []
        last_range = []
        for seq in ppb.build_peptides(struct):
            print seq.get_sequence()
        #read in the full sequence from the pdb file
        full_sequence = ''
        lines = fp.readlines()
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

        missing_set = set(range(first_missing[0], last_missing[0]))
        b = set(range(725, 727))
        c = list(missing_set & b)
        print c 


        