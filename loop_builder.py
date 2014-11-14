#!/usr/bin/env 
from modeller import *
from modeller.automodel import *
from Bio import PDB as pdb
import re
import csv
import os, glob

"""
loop_builder.py

This program fills in missing residues from (correctly formatted) pdb files.
It reads in the coordinates and detects any missing residues and then
gives an alignment file that can be used for loop modeling. It currently
does not do a very good job, as modeller is both poorly executed and also
I do not have the practical knowlede to make it do exactly what I want.

Known issues: 
    the paths are hard coded, you will need to modify them

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
        *** use the mutate_model.py script from modeller ***
    maybe put in a sanity check at the end that compares models of similar
        sequence or structure
"""

datafile = open('./structures.csv', 'r') #Opens the structures file for reading

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

#print pdb_info[0].id

for i in range(1, len(pdb_info)):
    pdb_name = pdb_info[i].id #saves given pdb name as a variable
    protein_name = pdb_info[i].protein #saves given protein name as a variable
    complete = pdb_info[i].complete #saves yes or no for complete
    structure_conf = pdb_info[i].conformation #saves active or inactive for conformation
    mutation = pdb_info[i].mutation
    print pdb_name
    pdb_file = './PDBs/'+pdb_name+'.pdb'
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
        #put in newlines into structure_sequence for proper PIR format
        for i in range(0,len(structure_sequence),70):
            structure_sequence = structure_sequence[:i] + "\n" + structure_sequence[i:]

        #read in the full sequence from the pdb file
        full_sequence = ''
        lines = fp.readlines()
        first = lines[0]
        header = re.split('HEADER\W+',first)[1] #uses a modified version of PDB file
        header_list = header.split(':')
        first_res = int(header_list[1])
        last_res = int(header_list[3])

        
        for index in range(1,10):
            split_line = re.split('REMARK 300 ',lines[index]) #appears that changing remark to 465 lowers the number of atoms
            if split_line[0] == '':
                full_sequence += split_line[1]


        #write the alignment file
        
        pdb_code = (pdb_name.split("-"[0])) 
        name = pdb_code[0] #changed this from hard coded 4F7S; it does not seem like this variable is used anywhere else
        PIR = open('active.ali','w')
        PIR.write(">P1;{0}\n".format(pdb_name))
        PIR.write("structureX:{0}".format(header))
        PIR.write("{0}*\n\n".format(structure_sequence.strip()))
        PIR.write(">P1;{0}\n".format(protein_name))
        PIR.write("sequence:{0}".format(header))
        PIR.write("{0}*\n\n".format(full_sequence.strip()))
        PIR.close()


        #begin modeller stuff here

        log.verbose()
        #where to look for pdb files
        
        env = environ()
        env.io.atom_files_directory = ['.', './PDBs']
        first_range = map(int, first_range) #makes this an integer array
        last_range = map(int, last_range) #makes this an integer array
        first_res_in_range = first_range.pop(0) #gets rid of the first element
        last_res_in_range = last_range.pop(-1) #gets rid of the last element
        first_missing = [x + 1 for x in last_range] #will use this to make missing residue ranges
        last_missing = [x - 1 for x in first_range] #will use this to make missing residue ranges

        
        # Create a new class based on 'loopmodel' so that we can redefine
        # select_loop_atoms
        class MyLoop(automodel):
            # This routine picks the residues to be refined by loop modeling
            def select_loop_atoms(self): #need to fix this
                #####make this easier to read by just doing if first case, add, if last case, add, etc.
                if last_res_in_range != last_res and first_res_in_range == first_res:
                    if not first_missing:
                        return selection(self.residue_range(last_res_in_range + ':', last_res + ':'))
                    else: 
                        return selection(self.residue_range(first_missing + ':', last_missing + ':'),
                                         self.residue_range(last_res_in_range + ':', last_res + ':'))
                if first_res_in_range != first_res and last_res_in_range == last_res:
                    if not last_missing:
                        return selection(self.residue_range(first_res + ':', first_res_in_range + ':'))
                    else: 
                        return selection(self.residue_range(first_res + ':', first_res_in_range + ':'),
                                         self.residue_range(first_missing + ':', last_missing + ':'))
                if first_res_in_range != first_res and last_res_in_range != last_res:
                    if not first_missing: # can use first_missing because if first_missing is empty, so is last_missing
                        return selection(self.residue_range(first_res + ':', first_res_in_range + ':'),
                                         self.residue_range(last_res_in_range + ':', last_res + ':'))
                    else:
                        return selection(self.residue_range(first_res + ':', first_res_in_range + ':'),
                                         self.residue_range(first_missing + ':', last_missing + ':'),
                                         self.residue_range(last_res_in_range + ':', last_res + ':'))
                else: #if first_res_in_range == first_res and last_res_in_range == last_res:
                    return selection(self.residue_range(first_missing + ':', last_missing + ':'))
                # index of the last model         
        
        a = MyLoop(env,
                alnfile  = 'active.ali',      # alignment filename
                knowns   = pdb_name,               # codes of the templates
                sequence = protein_name,               # code of the target
                library_schedule = autosched.slow,
                deviation = 1,
                assess_methods = assess.DOPE) # assess each loop with DOPE
        a.starting_model = 1                 # index of the first model
        a.ending_model  = 1 
        a.make() #do modeling and loop refinement

        #move directory and change name based on active/inactive and incomplete/complete designation
        if re.match("active", structure_conf) is not None:
            if re.match("yes", complete) is not None:
                os.rename(protein_name+'.B99990001.pdb', './actives/complete/'+protein_name+'_active') 
                os.rename(protein_name+'.D00000001', './actives/complete/'+protein_name+'_active_logFile')
            else:
                os.rename(protein_name+'.B99990001.pdb', './actives/incomplete/'+protein_name+'_active') 
                os.rename(protein_name+'.D00000001', './actives/incomplete/'+protein_name+'_active_logFile')
        if re.match("inactive", structure_conf) is not None: 
            if re.match("yes", complete) is not None:
                os.rename(protein_name+'.B99990001.pdb', './inactives/complete/'+protein_name+'_inactive') 
                os.rename(protein_name+'.D00000001', './inactives/complete/'+protein_name+'_inactive_logFile')
            else:
                os.rename(protein_name+'.B99990001.pdb', './inactives/incomplete/'+protein_name+'_inactive') 
                os.rename(protein_name+'.D00000001', './inactives/incomplete/'+protein_name+'_inactive_logFile')
        for filename in glob.glob("./"+protein_name+"*"):
            os.remove(filename) 