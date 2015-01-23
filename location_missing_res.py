#!/usr/bin/env 
from modeller import *
from modeller.optimizers import molecular_dynamics, conjugate_gradients
from modeller.automodel import *
from modeller.scripts import complete_pdb
from Bio import PDB as pdb
import re
import csv
import os, glob
import random
import numpy

"""
location_missing_res.py

This program uses all_kinase.csv in order to determine the location
of missing residues within kinase domain subdomains.  This is useful information, as
it will tell us the sections that loop_builder.py will be building back.  
"""

datafile = open('./structures2.csv', 'r') #Opens the structures file for reading
datareader = csv.reader(datafile) #reads structures file
data = [] #initializes a list called data
for row in datareader:
    data.append(row) #adds an element to data for each row in structures.csv 

datafile2 = open('./all_kinase.csv', 'r') #Opens the structures file for reading
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
        self.catstart = row[7]
        self.catend = row[8]
        self.activationloopstart = row[9]
        self.activationloopend = row[10]


pdb_info = [PDB_info(item) for item in data]
kinase_domain_info = [Kinase_domain_info(item) for item in data2] 
#print pdb_info[0].id

for i in range(1, len(pdb_info)):
    pdb_name = pdb_info[i].id #saves given pdb name as a variable
    protein_name = pdb_info[i].protein #saves given protein name as a variable
    complete = pdb_info[i].complete #saves yes or no for complete
    structure_conf = pdb_info[i].conformation #saves active or inactive for conformation
    mutation = pdb_info[i].mutation
    for j in range (1, len(kinase_domain_info)):
        if (protein_name == kinase_domain_info[j].protein_name):
            kd_start = int(kinase_domain_info[j].kdstart)
            kd_end = int(kinase_domain_info[j].kdend)
            ploop_start = int(kinase_domain_info[j].ploopstart)
            ploop_end = int(kinase_domain_info[j].ploopend)
            alphac_start = int(kinase_domain_info[j].alphacstart)
            alphac_end = int(kinase_domain_info[j].alphacend)
            cat_start = int(kinase_domain_info[j].catstart)
            cat_end = int(kinase_domain_info[j].catend)
            actloop_start = int(kinase_domain_info[j].activationloopstart)
            actloop_end = int(kinase_domain_info[j].activationloopend)
    in_kd = range(kd_start, kd_end + 1)
    in_ploop = range(ploop_start, ploop_end + 1)
    in_alphac = range(alphac_start, alphac_end + 1)
    in_cat = range(cat_start, cat_end + 1)
    in_actloop = range(actloop_start, actloop_end + 1)


    """
    search for a given string where protein_name == kinase_domain_info.protein_name
    then go to that line
    use regex to do this
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

        first_range = map(int, first_range) #makes this an integer array
        last_range = map(int, last_range) #makes this an integer array
        first_res_in_range = first_range.pop(0) #gets rid of the first element
        last_res_in_range = last_range.pop(-1) #gets rid of the last element
        first_missing = [x + 1 for x in last_range] #will use this to make missing residue ranges
        last_missing = [x - 1 for x in first_range] #will use this to make missing residue ranges

        ploop_counter = 0
        alphac_counter = 0
        cat_counter = 0
        actloop_counter = 0
        else_counter = 0

        if (first_res == kd_start and last_res == kd_end):
            print "good to go"
            for i in range(len(first_missing)):
                missing_set = range(first_missing[i], last_missing[i] + 1)
                for j in range(len(missing_set)):
                    if missing_set[j] in in_ploop:
                        ploop_counter += 1
                    if missing_set[j] in in_alphac:
                        alphac_counter += 1
                    if missing_set[j] in in_cat:
                        cat_counter += 1
                    if missing_set[j] in in_actloop:
                        actloop_counter += 1
                    else:
                        else_counter += 1
            missing_subdomains = open('./missing_subdomains.csv', 'a')
            missing_subdomains.write(protein_name+'_'+structure_conf+','+str(ploop_counter)+','+str(alphac_counter)+','+str(cat_counter)+','+str(actloop_counter)+','+str(else_counter)+"\n")
        else:
            print "looks like we've got a problem here"
            print first_res
            print last_res
            print kd_start
            print kd_end


        


        
