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
        
        if first_res != first_res_in_range:
            diff = first_res_in_range - first_res
            if (diff > 0):
                structure_sequence = diff*'-' + structure_sequence
        if last_res != last_res_in_range:
            diff = last_res - last_res_in_range
            if (diff > 0):
                structure_sequence = structure_sequence + diff*'-'

        #put in newlines into structure_sequence for proper PIR format
        for i in range(0,len(structure_sequence),70):
            structure_sequence = structure_sequence[:i] + "\n" + structure_sequence[i:]
        
        for index in range(1,10):
            split_line = re.split('REMARK 300 ',lines[index]) #appears that changing remark to 465 lowers the number of atoms
            if split_line[0] == '':
                full_sequence += split_line[1]


        #write the alignment file
        
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

        
        #begin modeller stuff here

        log.verbose()
        #where to look for pdb files
        
        env = environ()
        env.io.atom_files_directory = ['.', './PDBs']

        
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
                elif first_res_in_range != first_res and last_res_in_range == last_res:
                    if not last_missing:
                        return selection(self.residue_range(first_res + ':', first_res_in_range + ':'))
                    else: 
                        return selection(self.residue_range(first_res + ':', first_res_in_range + ':'),
                                         self.residue_range(first_missing + ':', last_missing + ':'))
                elif first_res_in_range != first_res and last_res_in_range != last_res:
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
        fp.close()
        #move directory and change name based on active/inactive and incomplete/complete designation
        if re.match("active", structure_conf) is not None:
            if re.match("yes", complete) is not None:
                os.rename(protein_name+'.B99990001.pdb', './actives/complete/'+protein_name+'_active.pdb') 
                os.rename(protein_name+'.D00000001', './actives/complete/'+protein_name+'_active_logFile')
                modelname = ('./actives/complete/'+protein_name+'_active.pdb')
                briefname = ('./actives/complete/'+protein_name+'_active')
                briefername = (protein_name+'_active')
                dir_name = ('./actives/complete/')
                is_active = 1
            else:
                os.rename(protein_name+'.B99990001.pdb', './actives/incomplete/'+protein_name+'_active.pdb') 
                os.rename(protein_name+'.D00000001', './actives/incomplete/'+protein_name+'_active_logFile')
                modelname = ('./actives/incomplete/'+protein_name+'_active.pdb')
                briefname = ('./actives/incomplete/'+protein_name+'_active')
                briefername = (protein_name+'_active')
                dir_name = ('./actives/incomplete/')
                is_active = 1
        if re.match("inactive", structure_conf) is not None: 
            if re.match("yes", complete) is not None:
                os.rename(protein_name+'.B99990001.pdb', './inactives/complete/'+protein_name+'_inactive.pdb') 
                os.rename(protein_name+'.D00000001', './inactives/complete/'+protein_name+'_inactive_logFile')
                modelname = ('./inactives/complete/'+protein_name+'_inactive.pdb')
                briefname = ('./inactives/complete/'+protein_name+'_inactive')
                briefername = (protein_name+'_inactive')
                dir_name = ('./inactives/complete/')
                is_active = 0
            else:
                os.rename(protein_name+'.B99990001.pdb', './inactives/incomplete/'+protein_name+'_inactive.pdb') 
                os.rename(protein_name+'.D00000001', './inactives/incomplete/'+protein_name+'_inactive_logFile')
                modelname = ('./inactives/incomplete/'+protein_name+'_inactive.pdb')
                briefname = ('./inactives/incomplete/'+protein_name+'_inactive')
                briefername = (protein_name+'_inactive')
                dir_name = ('./inactives/incomplete/')
                is_active = 0
        for filename in glob.glob("./"+protein_name+"*"):
            os.remove(filename)
         
        #mutate_model.py
     
        gromacs.editconf(f = modelname, resnr = first_res, o = modelname)
        

        if re.search('no_mutation', mutation) is not None:
            print 'No mutations here'
        else:
            print "Looks like we've got a mutation.  Let's check it out. \n"
            different_mutations = re.split('&', mutation)
            for mutant_res in range(0, len(different_mutations)):
                is_mutated = re.search(r"([a-z])([0-9]+)([a-z])", different_mutations[mutant_res], re.I)
                if is_mutated:
                    mutations_list = is_mutated.groups()
                    respos = mutations_list[1]
                    restyp = pdb.Polypeptide.one_to_three(mutations_list[0]) #get three letter code

                    #makes use of the optimize function in modeller
                    def optimize(atmsel, sched):
                        #conjugate gradient
                        for step in sched:
                            step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
                        #md
                        refine(atmsel)
                        cg = conjugate_gradients()
                        cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


                    #molecular dynamics
                    def refine(atmsel):
                        # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
                        md = molecular_dynamics(cap_atom_shift=0.39, md_time_step=4.0,
                                                md_return='FINAL')
                        init_vel = True
                        for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                                    (200, 600,
                                                     (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
                            for temp in temps:
                                md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                                             max_iterations=its, equilibrate=equil)
                                init_vel = False


                    #use homologs and dihedral library for dihedral angle restraints
                    def make_restraints(mdl1, aln):
                       rsr = mdl1.restraints
                       rsr.clear()
                       s = selection(mdl1)
                       for typ in ('stereo', 'phi-psi_binormal'):
                           rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
                       for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
                           rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                                    spline_dx=0.3, spline_min_points = 5, aln=aln,
                                    spline_on_site=True)

                    log.verbose()

                    # Set a different value for rand_seed to get a different final model
                    env = environ(rand_seed=-49837)

                    env.io.hetatm = True
                    #soft sphere potential
                    env.edat.dynamic_sphere=False
                    #lennard-jones potential (more accurate)
                    env.edat.dynamic_lennard=True
                    #https://salilab.org/modeller/manual/node127.html to learn more about contact_shell and update_dynamic
                    env.edat.contact_shell = 4.0
                    env.edat.update_dynamic = 0.39

                    # Read customized topology file with phosphoserines (or standard one)
                    env.libs.topology.read(file='$(LIB)/top_heav.lib')

                    # Read customized CHARMM parameter library with phosphoserines (or standard one)
                    env.libs.parameters.read(file='$(LIB)/par.lib')


                    # Read the original PDB file and copy its sequence to the alignment array:
                    mdl1 = model(env, file=modelname)
                    ali = alignment(env)
                    ali.append_model(mdl1, atom_files=modelname, align_codes=modelname)

                    #set up the mutate residue selection segment
                    s = selection(mdl1.residues[respos])

                    #perform the mutate residue operation
                    s.mutate(residue_type=restyp)
                    #get two copies of the sequence.  A modeller trick to get things set up
                    ali.append_model(mdl1, align_codes=modelname)

                    # Generate molecular topology for mutant
                    mdl1.clear_topology()
                    mdl1.generate_topology(ali[-1])


                    # Transfer all the coordinates you can from the template native structure
                    # to the mutant (this works even if the order of atoms in the native PDB
                    # file is not standard):
                    #here we are generating the model by reading the template coordinates
                    mdl1.transfer_xyz(ali)

                    # Build the remaining unknown coordinates
                    mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

                    #yes model2 is the same file as model1.  It's a modeller trick.
                    mdl2 = model(env, file=modelname)

                    #required to do a transfer_res_numb
                    #ali.append_model(mdl2, atom_files=modelname, align_codes=modelname)
                    #transfers from "model 2" to "model 1"
                    mdl1.res_num_from(mdl2,ali)

                    #It is usually necessary to write the mutated sequence out and read it in
                    #before proceeding, because not all sequence related information about MODEL
                    #is changed by this command (e.g., internal coordinates, charges, and atom
                    #types and radii are not updated).

                    mdl1.write(file=modelname+restyp+respos+'.tmp')
                    mdl1.read(file=modelname+restyp+respos+'.tmp')

                    #set up restraints before computing energy
                    #we do this a second time because the model has been written out and read in,
                    #clearing the previously set restraints
                    make_restraints(mdl1, ali)

                    #a non-bonded pair has to have at least as many selected atoms
                    mdl1.env.edat.nonbonded_sel_atoms=1

                    sched = autosched.loop.make_for_model(mdl1)

                    #only optimize the selected residue (in first pass, just atoms in selected
                    #residue, in second pass, include nonbonded neighboring atoms)
                    #set up the mutate residue selection segment
                    s = selection(mdl1.residues[respos])

                    mdl1.restraints.unpick_all()
                    mdl1.restraints.pick(s)

                    s.energy()

                    s.randomize_xyz(deviation=4.0)

                    mdl1.env.edat.nonbonded_sel_atoms=2
                    optimize(s, sched)

                    #feels environment (energy computed on pairs that have at least one member
                    #in the selected)
                    mdl1.env.edat.nonbonded_sel_atoms=1
                    optimize(s, sched)

                    s.energy()

                    #give a proper name
                    mdl1.write(file=briefname+restyp+respos+'.pdb')

                    #delete the temporary file
                    os.remove(modelname+restyp+respos+'.tmp')
                    
                    modelname = (briefname+restyp+respos+'.pdb')
                    briefname = (briefname+restyp+respos)
                    briefername = (briefername+restyp+respos)
                else:
                    print 'Mutation of '+pdb_name+' not recognized. Going to take note of this.'
                    bad_mutation = open('bad_mutation_PDBs.csv','a')
                    bad_mutation.write(pdb_name+','+protein_name+','+complete+','+structure_conf+','+mutation+"\n")
                    bad_mutation.close()
        if (is_active == 1):
            template = ('./actives/complete/ABL1_active.pdb')
            template_briefname = ('ABL1_active')
            pml_viewer = open('./active_aligner.pml', 'a')
            pml_viewer.write('load '+modelname+"\n"+'hide lines, '+briefername+"\n"+'show cartoon, '
                             +briefername+"\n"+'align '+briefername+', '+template_briefname+', '+'cycles=0'"\n"+"\n")
        else:
            template = ('./inactives/complete/ABL1_inactive.pdb')
            template_briefname = ('ABL1_inactive')
            pml_viewer = open('./inactive_aligner.pml', 'a')
            pml_viewer.write('load '+modelname+"\n"+'hide lines, '+briefername+"\n"+'show cartoon, '
                             +briefername+"\n"+'align '+briefername+', '+template_briefname+', '+'cycles=0'"\n"+"\n")
   
#write the align-multiple.ali
        """
        env = environ()
        env.io.atom_files_directory = ['./']
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        
        template = ('./actives/complete/ABL1_active.pdb')
        mdl = complete_pdb(env, modelname)
        mdl2 = complete_pdb(env, template)
        # Write out chain sequences:
        for c in mdl.chains:
            c.write(file='align-multiple1.ali', atom_file=modelname,
                    align_code=modelname)
        for c in mdl2.chains:
            c.write(file='align-multiple2.ali', atom_file=template,
                    align_code=template)
        filenames = ['align-multiple1.ali', 'align-multiple2.ali']
        with open('./align-multiple.ali', 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    outfile.write(infile.read())

        log.verbose()
        env = environ()
        env.io.atom_files_directory = ['./']

        aln = alignment(env)
        for (code) in ((template), ('1uld', 'D'), ('1ulf', 'B'),
                              ('1ulg', 'B'), ('1is5', 'A')):
            mdl = model(env, file=code, model_segment=('FIRST:A', 'LAST:A'))
            aln.append_model(mdl, atom_files=code, align_codes=code)

        for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                            ((1., 0.5, 1., 1., 1., 0.), False, True),
                                            ((1., 1., 1., 1., 1., 0.), True, False)):
            aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                       rr_file='$(LIB)/as1.sim.mat', overhang=30,
                       gap_penalties_1d=(-450, -50),
                       gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
                       dendrogram_file='1is3A.tree',
                       alignment_type='tree', # If 'progresive', the tree is not
                                              # computed and all structues will be
                                              # aligned sequentially to the first
                       #ext_tree_file='1is3A_exmat.mtx', # Tree building can be avoided
                                                         # if the tree is input
                       feature_weights=weights, # For a multiple sequence alignment only
                                                # the first feature needs to be non-zero
                       improve_alignment=True, fit=True, write_fit=write_fit,
                       write_whole_pdb=whole, output='ALIGNMENT QUALITY')

        aln.write(file='1is3A.pap', alignment_format='PAP')
        aln.write(file='1is3A.ali', alignment_format='PIR')

        # The number of equivalent positions at different RMS_CUTOFF values can be
        # computed by changing the RMS value and keeping all feature weights = 0
        aln.salign(rms_cutoff=1.0,
                   normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30,
                   gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
                   gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
                   alignment_type='progressive', feature_weights=[0]*6,
                   improve_alignment=False, fit=False, write_fit=True,
                   write_whole_pdb=False, output='QUALITY')
        

        log.verbose()    # request verbose output
        env = environ()  # create a new MODELLER environment to build this model in

        # directories for input atom files
        env.io.atom_files_directory = ['./']

        a = automodel(env,
                      alnfile  = 'align-multiple.ali', # alignment filename
                      knowns   = (template),     # codes of the templates
                      sequence = modelname)               # code of the target
        a.starting_model= 1                 # index of the first model
        a.ending_model  = 1                 # index of the last model
                                            # (determines how many models to calculate)
        a.make()                            # do the actual comparative modeling
        """