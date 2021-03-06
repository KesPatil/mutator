import math
import scipy
import scipy.constants
import scipy.stats as spy
import numpy as np
import matplotlib.pyplot as plt
import random
import os
from heapq import *

POS_X_START = 31-1
POS_X_END = 38-1
POS_Y_START = 39-1
POS_Y_END = 46-1
POS_Z_START = 47-1
POS_Z_END = 54-1
RESIDUE_NUM_START = 23-1
RESIDUE_NUM_END = 26-1
ACTIVE_START = 17
INACTIVE_START = 19
ACTIVE_END = 11
INACTIVE_END = 13
ACTIVE_START_INCOMP = 19
INACTIVE_START_INCOMP = 21
AMINO_START = 18-1
AMINO_END = 20-1
ATOM_START = 13-1
ATOM_END = 16-1

# Functions for finding the array of minimums, array of maximums, and array of means of all rows in a 2-D array
def min_array(array_2d):
    mins = []
    for a in array_2d:
        mins.append(min(a))
    return mins

def max_array(array_2d):
    maxs = []
    for a in array_2d:
        maxs.append(max(a))
    return maxs

def mean_array(array_2d):
    means = []
    for a in array_2d:
        means.append(sum(a)/len(a))
    return means

# Give the name of a file containing information about kinases (e.g. residue numbers for the activation loop, catalytic loop, nucleotide-binding loop, and alpha-c helix), return a dictionary with the name of the kinase as the key and the values of interest
def extract_endpoints(file_name):
    start_end_file = open(file_name, "r+")
    start_end_lines = start_end_file.readlines()
    endpoints = {}
    for l in range(1, len(start_end_lines)):
        line_list = start_end_lines[l].split(",")
        endpoints[line_list[0]] = []
        for l2 in range(1, len(line_list)):
            endpoints[line_list[0]].append(line_list[l2])
    return endpoints

# Divide all available data (contained within two files with file names ac_fn and inac_fn) randomly into a training and a testing set, where the proportion of the data (where the size of the data is measured by the smaller of the active and inactive datasets) that is put into the training set is set by a parameter within the function. The flag variable indicates whether all the data should be put into the training set (flag=1) or not (flag=0).
def get_training_testing_sets(ac_fn, inac_fn, flag):
    ac_file = open(ac_fn, "r+")
    inac_file = open(inac_fn, "r+")
    ac_lines = ac_file.readlines()
    inac_lines = inac_file.readlines()
    training = []
    testing = []
    proportion = 0.5
    if (flag == 1):
        for acl in ac_lines:
            training.append([1, acl[:len(acl)-1]])
        for inacl in inac_lines:
            training.append([0, inacl[:len(inacl)-1]])
    else:
        train_count = int(min(len(ac_lines), len(inac_lines))*proportion)
        act_trains = []
        inact_trains = []
        for i in range(train_count):
            num = random.randint(0, len(ac_lines)-1)
            while (num in act_trains):
                num = random.randint(0, len(ac_lines)-1)
            act_trains.append(num)
            num = random.randint(0, len(inac_lines)-1)
            while (num in inact_trains):
                num = random.randint(0, len(inac_lines)-1)
            inact_trains.append(num)
        for j in range(0, len(ac_lines)):
            if (j in act_trains):
                training.append([1, ac_lines[j][:len(ac_lines[j])-1]])
            else:
                testing.append([1, ac_lines[j][:len(ac_lines[j])-1]])
        for k in range(0, len(inac_lines)):
            if (k in inact_trains):
                training.append([0, inac_lines[k][:len(inac_lines[k])-1]])
            else:
                testing.append([0, inac_lines[k][:len(inac_lines[k])-1]])
    print str(len(ac_lines)) + ", " + str(len(inac_lines))
    print str(len(training)) + ", " + str(len(testing))
    return [training, testing]

# Tokenize a line in a pdb file
def pdb_tokenize(pdb_line):
    words =[]
    curr_word = pdb_line[0]
    for i in range(1, len(pdb_line)):
        if (pdb_line[i] != " "):
            curr_word += pdb_line[i]
        elif (len(curr_word) >= 1):
            words.append(curr_word)
            curr_word = ""
    if (len(curr_word) >= 2):
        words.append(curr_word[:len(curr_word)-1])
    return words

# Given two lists/tuples/etc. that represent two points in 3-space, return the Euclidean distance between the two points
def dist_3d(xyz1, xyz2):
    return math.sqrt((xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 + (xyz1[2]-xyz2[2])**2)

# Given a list of the lines in a pdb (pdb_lines), the index of the line in the pdb that corresponds to a particular atom (single), the starting and ending residue numbers for a specific region in the kinase such as the catalytic loop (list_start and list_end), and the index of the line in the pdb to start searching for the residues in that region (index)
def dist_list(pdb_lines, single, index, list_start, list_end):
    single_x = float(pdb_lines[single][POS_X_START:POS_X_END+1])
    single_y = float(pdb_lines[single][POS_Y_START:POS_Y_END+1])
    single_z = float(pdb_lines[single][POS_Z_START:POS_Z_END+1])
    #print str(single) + " " + str(single_x) + " " + str(single_y) + " " + str(single_z)
    coords1 = [single_x, single_y, single_z]
    d_list = []
    j = index
    while (j < len(pdb_lines) and "ATOM" in pdb_lines[j]):
        #print pdb_lines[j][RESIDUE_NUM_START:RESIDUE_NUM_END+1]
        if (int(pdb_lines[j][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= list_start and int(pdb_lines[j][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= list_end):
            this_elements = pdb_lines[j]
            coords2 = [float(this_elements[POS_X_START:POS_X_END+1]), float(this_elements[POS_Y_START:POS_Y_END+1]), float(this_elements[POS_Z_START:POS_Z_END+1])]
            d_list.append(dist_3d(coords1, coords2))
        j += 1
    return d_list

# Find all pairwise distances between 'CA' atoms in two different regions such as the activation and catalytic loops. Given the name of the kinase (kinase_name), a list of the lines in the pdb (pdb_lines), a list/tuple of the starting and ending residue numbers for the first region (region1), and a list/tuple of the starting and ending residue numbers for the second region (region2), return a 2-D array.
def two_region_distances(kinase_name, pdb_lines, region1, region2):
    reg1_start = region1[0]
    reg1_end = region1[1]
    reg2_start = region2[0]
    reg2_end = region2[1]
    i = 0
    ind = -1
    distances = []
    while (i < len(pdb_lines)):
        if ("ATOM" in pdb_lines[i]):
            if (ind == -1):
                ind = i
            if (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= reg1_start and int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= reg1_end):
                if (pdb_lines[i][ATOM_START:ATOM_END+1].strip(' ') == "CA"):
                    the_dist_list = dist_list(pdb_lines, i, ind, reg2_start, reg2_end)
                    distances.append(the_dist_list)
        i += 1
    return distances

# Compute the minimum distance between any atom in the activation loop and any atom in the catalytic loop. This function returns the correlation between that quantity (in addition to correlations for the max-minimum distance and the mean-minimum distance) and whether a kinase is active and inactive given a list of kinases (training). Should use two_region_distances() with min_array() instead of this function.
def min_act_cat_dist(training, endpoints_file):
    min_active_list = []
    min_inactive_list = []
    max_active_list = []
    max_inactive_list = []
    mean_active_list = []
    mean_inactive_list = []
    endpoint_info = extract_endpoints(endpoints_file)
    for t in training:
        pdb_file = open(t[1], "r+")
        pdb_lines = pdb_file.readlines()
        kinase_name = ""
        if (t[0] == 1):
            kinase_name = t[1][ACTIVE_START:len(t[1])-ACTIVE_END]
        else:
            kinase_name = t[1][INACTIVE_START:len(t[1])-INACTIVE_END]
        if (kinase_name not in endpoint_info):
            continue
        kinase_endpoints = endpoint_info[kinase_name]
        cat_start = int(kinase_endpoints[len(kinase_endpoints)-4])
        cat_end = int(kinase_endpoints[len(kinase_endpoints)-3])
        act_start = int(kinase_endpoints[len(kinase_endpoints)-2])
        act_end = int(kinase_endpoints[len(kinase_endpoints)-1])
        flag = 0
        ind = 0
        i = 0
        min_distances = []
        max_distances = []
        mean_distances = []
        while (i < len(pdb_lines)):
            if ("ATOM" in pdb_lines[i]):
                if (flag == 0):
                    flag = 1
                    ind = i
                if (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= act_start and int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= act_end):
                    #print t[1]
                    the_dist_list = dist_list(pdb_lines, i, ind, cat_start, cat_end)
                    min_dist_list = min(the_dist_list)
                    max_dist_list = max(the_dist_list)
                    mean_dist_list = sum(the_dist_list)/len(the_dist_list)
                    min_distances.append(min_dist_list)
                    max_distances.append(max_dist_list)
                    mean_distances.append(mean_dist_list)
            i += 1
        if (t[0] == 1):
            min_active_list.append(sum(min_distances)/len(min_distances))
            max_active_list.append(sum(max_distances)/len(max_distances))
            mean_active_list.append(sum(mean_distances)/len(mean_distances))
        else:
            min_inactive_list.append(sum(min_distances)/len(min_distances))
            max_inactive_list.append(sum(max_distances)/len(max_distances))
            mean_inactive_list.append(sum(mean_distances)/len(mean_distances))
    classifications = []
    min_datasets = []
    max_datasets = []
    mean_datasets = []
    for a in range(0, len(min_active_list)):
        classifications.append(1)
        min_datasets.append(min_active_list[a])
        max_datasets.append(max_active_list[a])
        mean_datasets.append(mean_active_list[a])
    for ina in range(0, len(min_inactive_list)):
        classifications.append(0)
        min_datasets.append(min_inactive_list[ina])
        max_datasets.append(max_inactive_list[ina])
        mean_datasets.append(mean_inactive_list[ina])
    #print "Using min distances:"
    #print spy.pearsonr(classifications, min_datasets)
    #print "Using max distances:"
    #print spy.pearsonr(classifications, max_datasets)
    #print "Using mean distances:"
    #print spy.pearsonr(classifications, mean_datasets)
    #return [active_list, inactive_list]
    return [spy.pearsonr(classifications, min_datasets)[0], spy.pearsonr(classifications, max_datasets)[0], spy.pearsonr(classifications, mean_datasets)[0]]

# Given a list of lines in a pdb and an index of a line in the pdb, return a pair (i.e. a list of size 2) where the first element is the name of the next amino acid after the given index and the second element is the index where this next amino acid begins.
def find_next_amino(pdb_lines, currIndex):
    if ("ATOM" not in pdb_lines[currIndex]):
        return ["", -1]
    currAmino = pdb_lines[currIndex][AMINO_START:AMINO_END+1]
    while (pdb_lines[currIndex][AMINO_START:AMINO_END+1] == currAmino):
        currIndex += 1
    return [pdb_lines[currIndex][AMINO_START:AMINO_END+1], currIndex]

# Compute the distance between the atoms in the DFG residue (in the activation loop) and the atoms in the catalytic loop for several active and inactive residues. Return three correlations between whether or not a kinase is active and 1) minimum distance, 2) maximum distance, and 3) mean distance.
def avg_dfg_dist(training, endpoints_file):
    min_active_list = []
    max_active_list = []
    mean_active_list = []
    min_inactive_list = []
    max_inactive_list = []
    mean_inactive_list = []
    endpoint_info = extract_endpoints(endpoints_file)
    for t in training:
        print t[1]
        pdb_file = open(t[1], "r+")
        pdb_lines = pdb_file.readlines()
        kinase_name = ""
        if (t[0] == 1):
            kinase_name = t[1][ACTIVE_START:len(t[1])-ACTIVE_END]
        else:
            kinase_name = t[1][INACTIVE_START:len(t[1])-INACTIVE_END]
        if (kinase_name not in endpoint_info):
            continue
        kinase_endpoints = endpoint_info[kinase_name]
        cat_start = int(kinase_endpoints[len(kinase_endpoints)-4])
        cat_end = int(kinase_endpoints[len(kinase_endpoints)-3])
        act_start = int(kinase_endpoints[len(kinase_endpoints)-2])
        act_end = int(kinase_endpoints[len(kinase_endpoints)-1])
        i = 0
        ind = 0
        flag = 0
        min_distances = []
        max_distances = []
        mean_distances = []
        while (i < len(pdb_lines)):
            if ("ATOM" in pdb_lines[i]):
                if (flag == 0):
                    flag = 1
                    ind = i
                if(int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= act_start and int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= act_end):
                    print pdb_lines[i][AMINO_START:AMINO_END+1]
                    if (pdb_lines[i][AMINO_START:AMINO_END+1] == "ASP"):
                        print "ASP" + str(i)
                        [i2_amino, i2] = find_next_amino(pdb_lines, i)
                        if (i2_amino == "PHE"):
                            print "PHE" + str(i)
                            [i3_amino, i3] = find_next_amino(pdb_lines, i2)
                            if (i3_amino == "GLY"):
                                [i4_amino, i4] = find_next_amino(pdb_lines, i3)
                                for j in range(i, i4):
                                    the_dist_list = dist_list(pdb_lines, j, ind, cat_start, cat_end)
                                    min_dist_list = min(the_dist_list)
                                    max_dist_list = max(the_dist_list)
                                    mean_dist_list = float(sum(the_dist_list))/len(the_dist_list)
                                    min_distances.append(min_dist_list)
                                    max_distances.append(max_dist_list)
                                    mean_distances.append(mean_dist_list)
                                break
            i += 1
        if (len(min_distances) > 0):
            if (t[0] == 1):
                min_active_list.append(sum(min_distances)/len(min_distances))
                max_active_list.append(sum(max_distances)/len(max_distances))
                mean_active_list.append(sum(mean_distances)/len(mean_distances))
            else:
                min_inactive_list.append(sum(min_distances)/len(min_distances))
                max_inactive_list.append(sum(max_distances)/len(max_distances))
                mean_inactive_list.append(sum(mean_distances)/len(mean_distances))
    classifications = []
    min_datasets = []
    max_datasets = []
    mean_datasets = []
    for a in range(0, len(min_active_list)):
        classifications.append(1)
        min_datasets.append(min_active_list[a])
        max_datasets.append(max_active_list[a])
        mean_datasets.append(mean_active_list[a])
    for ina in range(0, len(min_inactive_list)):
        classifications.append(0)
        min_datasets.append(min_inactive_list[ina])
        max_datasets.append(max_inactive_list[ina])
        mean_datasets.append(mean_inactive_list[ina])
    return [spy.pearsonr(classifications, min_datasets)[0], spy.pearsonr(classifications, max_datasets)[0], spy.pearsonr(classifications, mean_datasets)[0]]

# A crude measure of the degree to which a set of points in 3-space can be modeled by a plane. Computes the least-squares plane for the points and the residuals for the fit. Returns the correlation between the mean residual and whether a kinase is active or inactive
def planar_similarity(training, endpoints_file):
    active_list = []
    inactive_list = []
    endpoint_info = extract_endpoints(endpoints_file)
    dataset = []
    classifications = []
    for t in training:
        print t[1]
        pdb_file = open(t[1], "r+")
        pdb_lines = pdb_file.readlines()
        kinase_name = ""
        if (t[0] == 1):
            kinase_name = t[1][ACTIVE_START:len(t[1])-ACTIVE_END]
        else:
            kinase_name = t[1][INACTIVE_START:len(t[1])-INACTIVE_END]
        if (kinase_name not in endpoint_info):
            print "kinase name not in endpoint info"
            continue
        kinase_endpoints = endpoint_info[kinase_name]
        #act_start = int(kinase_endpoints[len(kinase_endpoints)-2])
        #act_end = int(kinase_endpoints[len(kinase_endpoints)-1])
        alphac_start = int(kinase_endpoints[len(kinase_endpoints)-6])
        alphac_end = int(kinase_endpoints[len(kinase_endpoints)-5])
        i = 0
        flag = 0
        X = []
        Z = []
        while (i < len(pdb_lines)):
            if ("ATOM" in pdb_lines[i]):
                #if (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= act_start) and (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= act_end):
                if (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= alphac_start) and (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= alphac_end):
                    x_value = float(pdb_lines[i][POS_X_START:POS_X_END+1])
                    y_value = float(pdb_lines[i][POS_Y_START:POS_Y_END+1])
                    z_value = float(pdb_lines[i][POS_Z_START:POS_Z_END+1])
                    X.append([x_value, y_value, 1])
                    Z.append(z_value)
                    flag = 1
                elif flag == 1:
                    break
            i += 1
        X = np.array(X)
        Z = np.array(Z)
        results = np.linalg.lstsq(X, Z)
        dataset.append(np.mean(results[1]))
        classifications.append(t[0])
    return spy.pearsonr(dataset, classifications)[0]

# Compute the electric field at a position due to a charge given the positions of the charge and the position of interest and the value of the charge by Coulomb's Law.
def electric_field(pos1, charge2, pos2):
    r_vec = pos2 - pos1
    mult_constant = 1.0/(4*scipy.constants.pi*scipy.constants.epsilon_0)
    return mult_constant*charge2*r_vec/(np.linalg.norm(r_vec)**3)

# Compute the electric potential at a position due to a charge given the positions of the charge and the location of interest and the value of the charge by Coulomb's Law
def electric_potential(pos1, charge2, pos2):
    r_vec = pos2-pos1
    mult_constant = 1.0/(4*scipy.constants.pi*scipy.constants.epsilon_0)
    return mult_constant*charge2/np.linalg.norm(r_vec)

# Compute correlations between whether a kinase is active or inactive and a rough estimate of the contribution to the electric field at the center of the catalytic loop due to the atoms in the activation loop.
def act_cat_field(training, endpoints_file):
    dataset_field = []
    dataset_potential = []
    classifications = []
    endpoint_info = extract_endpoints(endpoints_file)
    for t in training:
        print t[1]
        pdb_file = open(t[1], "r+")
        pdb_lines = pdb_file.readlines()
        kinase_name = ""
        if (t[0] == 1):
            kinase_name = t[1][ACTIVE_START:len(t[1])-ACTIVE_END]
        else:
            kinase_name = t[1][INACTIVE_START:len(t[1])-INACTIVE_END]
        if (kinase_name not in endpoint_info):
            print kinase_name + " not in endpoint_info"
            continue
        kinase_endpoints = endpoint_info[kinase_name]
        cat_start = int(kinase_endpoints[len(kinase_endpoints)-4])
        cat_end = int(kinase_endpoints[len(kinase_endpoints)-3])
        act_start = int(kinase_endpoints[len(kinase_endpoints)-2])
        act_end = int(kinase_endpoints[len(kinase_endpoints)-1])
        cat_center_pos = [0.0, 0.0, 0.0]
        i = 0
        flag = 0
        while (i < len(pdb_lines)):
            if ("ATOM" in pdb_lines[i]):
                if ((int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= cat_start) and (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= cat_end)):
                    flag = 1
                    x_value = float(pdb_lines[i][POS_X_START:POS_X_END+1])
                    y_value = float(pdb_lines[i][POS_Y_START:POS_Y_END+1])
                    z_value = float(pdb_lines[i][POS_Z_START:POS_Z_END+1])
                    cat_center_pos = [cat_center_pos[0]+x_value, cat_center_pos[1]+y_value, cat_center_pos[2]+z_value]
                elif (flag == 1):
                    break
            i += 1
        cat_center_pos = np.array(cat_center_pos)/(cat_end-cat_start+1)
        i = 0
        flag = 0
        cat_field = np.array([0.0, 0.0, 0.0])
        cat_potential = 0.0
        while (i < len(pdb_lines)):
            if ("ATOM" in pdb_lines[i]):
                if ((int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= act_start) and (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= act_end)):
                    flag = 1
                    x_value = float(pdb_lines[i][POS_X_START:POS_X_END+1])
                    y_value = float(pdb_lines[i][POS_Y_START:POS_Y_END+1])
                    z_value = float(pdb_lines[i][POS_Z_START:POS_Z_END+1])
                    act_pos = np.array([x_value, y_value, z_value])
                    charge = 0.0
                    atom_name = pdb_lines[i][ATOM_START:ATOM_END+1]
                    if ("NH1" in atom_name or "NH2" in atom_name or "NZ" in atom_name):
                        charge = 1.0
                    elif ("OE1" in atom_name or "OE2" in atom_name or "OD1" in atom_name or "OD2" in atom_name):
                        charge = -1.0
                    cat_field += electric_field(cat_center_pos, charge, act_pos)
                    cat_potential += electric_potential(cat_center_pos, charge, act_pos)
                elif (flag == 1):
                    break
            i += 1
        dataset_field.append(np.linalg.norm(cat_field))
        dataset_potential.append(cat_potential)
        classifications.append(t[0])
    return [spy.pearsonr(dataset_field, classifications)[0], spy.pearsonr(dataset_potential, classifications)[0]]

# Find the average correlation computed by min_act_cat_dist() over n trials.
def averageCorrelations():
    minSum = 0
    maxSum = 0
    meanSum = 0
    n = 100
    for i in range(n):
        [training, testing] =  get_training_testing_sets("activePDBs.txt", "inactivePDBs.txt")
        arr = min_act_cat_dist(training, "all_kinase.csv")
        minSum += arr[0]
        maxSum += arr[1]
        meanSum += arr[2]
    minSum = float(minSum)/float(n)
    maxSum = float(maxSum)/float(n)
    meanSum = float(meanSum)/float(n)
    return [minSum, maxSum, meanSum]

# For a given position (pos) in a given kinase (kinase_name), the function returns the align number associated with that residue.
def find_align_num(align_file, kinase_name, pos):
    align_f = open(align_file, "r+")
    align_info = align_f.readlines()
    align_num = -1
    running_align_num = -1
    running_align_line = -1
    for j in range(0, len(align_info)):
        split_line = align_info[j].split(" ")
        if (len(split_line) == 3):
            if (split_line[0] == "align"):
                running_align_num = int(split_line[2][:len(split_line[2])-1])
                running_align_line = j
            else:
                if (split_line[0][1:] == kinase_name and int(split_line[2][:len(split_line[2])-1]) == pos):
                    align_num = int(split_line[2][:len(split_line[2])-1])
                    break
    return align_num, running_align_line

def findSpecificDistance(file_name, endpoints_file, align_file, kinase_name, pos):
    list_file = open(file_name, "r+")
    list_lines = list_file.readlines()
    kinases = []
    endpoint_info = extract_endpoints(endpoints_file)
    i = 0
    while (i < len(list_lines)-1):
        line1 = list_lines[i]
        line2 = list_lines[i+1]
        line1_elements = line1.split(",")
        line2_elements = line2.split(",")
        if (line1_elements[1] == line2_elements[1]):
            to_add = []
            if (line1_elements[2] == "yes"):
                to_add.append("actives/complete/" + line1_elements[1] + "_active.pdb")
            else:
                to_add.append("actives/incomplete/" + line1_elements[1] + "_active.pdb")
            if (line2_elements[2] == "yes"):
                to_add.append("inactives/complete/" + line2_elements[1] + "_inactive.pdb")
            else:
                to_add.append("inactives/incomplete/" + line2_elements[1] + "_inactive.pdb")
            kinases.append(to_add)
            i += 2
        else:
            i += 1
    differences = []
    align_f = open(align_file, "r+")
    align_info = align_f.readlines()
    align_num, running_align_line = find_align_num(align_file, kinase_name, pos)
    active_dists = []
    inactive_dists = []
    for k in kinases:
        if (not os.path.isfile(k[0])) or (not os.path.isfile(k[1])): 
            continue
        active_file = open(k[0], "r+")
        inactive_file = open(k[1], "r+")
        active_lines = active_file.readlines()
        inactive_lines = inactive_file.readlines()
        k_name = ""
        if ("incomplete" in k[0]):
            k_name = k[0][ACTIVE_START_INCOMP:len(k[0])-ACTIVE_END]
        else:
            k_name = k[0][ACTIVE_START:len(k[0])-ACTIVE_END]
        kinase_endpoints = endpoint_info[k_name]
        cat_start = int(kinase_endpoints[len(kinase_endpoints)-4])
        cat_end = int(kinase_endpoints[len(kinase_endpoints)-3])
        i = 0
        active_dist = 0
        inactive_dist = 0
        equivalent_pos = -1
        act_mol_center = [0.0, 0.0, 0.0]
        cat_center = [0.0, 0.0, 0.0]
        count1 = 0
        count2 = 0
        for j in range(running_align_line, len(align_info)):
            split_line = align_info[j].split(" ")
            if (len(split_line) == 3):
                if (split_line[0] != "align"):
                    if (split_line[0][1:] == k_name):
                        equivalent_pos = int(split_line[2][:len(split_line[2])-1])
                        break
        while (i < len(active_lines)):
            if ("ATOM" in active_lines[i]):
                if (int(active_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1] ) == equivalent_pos):
                    x_value = float(active_lines[i][POS_X_START:POS_X_END+1])
                    y_value = float(active_lines[i][POS_Y_START:POS_Y_END+1])
                    z_value = float(active_lines[i][POS_Z_START:POS_Z_END+1])
                    act_mol_center[0] += x_value
                    act_mol_center[1] += y_value
                    act_mol_center[2] += z_value
                    count1 += 1
                elif (int(active_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= cat_start and int(active_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= cat_end):
                    x_value = float(active_lines[i][POS_X_START:POS_X_END+1])
                    y_value = float(active_lines[i][POS_Y_START:POS_Y_END+1])
                    z_value = float(active_lines[i][POS_Z_START:POS_Z_END+1])
                    cat_center[0] += x_value
                    cat_center[1] += y_value
                    cat_center[2] += z_value
                    count2 += 1
            i += 1
        act_mol_center[0] /= count1
        act_mol_center[1] /= count1
        act_mol_center[2] /= count1
        cat_center[0] /= count2
        cat_center[1] /= count2
        cat_center[2] /= count2
        active_dist = dist_3d(act_mol_center, cat_center)
        i = 0
        act_mol_center = [0.0, 0.0, 0.0]
        cat_center = [0.0, 0.0, 0.0]
        count1 = 0
        count2 = 0
        while (i < len(inactive_lines)):
            if ("ATOM" in inactive_lines[i]):
                if (int(inactive_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1] ) == equivalent_pos):
                    x_value = float(inactive_lines[i][POS_X_START:POS_X_END+1])
                    y_value = float(inactive_lines[i][POS_Y_START:POS_Y_END+1])
                    z_value = float(inactive_lines[i][POS_Z_START:POS_Z_END+1])
                    act_mol_center[0] += x_value
                    act_mol_center[1] += y_value
                    act_mol_center[2] += z_value
                    count1 += 1
                elif (int(inactive_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= cat_start and int(inactive_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= cat_end):
                    x_value = float(inactive_lines[i][POS_X_START:POS_X_END+1])
                    y_value = float(inactive_lines[i][POS_Y_START:POS_Y_END+1])
                    z_value = float(inactive_lines[i][POS_Z_START:POS_Z_END+1])
                    cat_center[0] += x_value
                    cat_center[1] += y_value
                    cat_center[2] += z_value
                    count2 += 1
            i += 1

        act_mol_center[0] /= count1
        act_mol_center[1] /= count1
        act_mol_center[2] /= count1
        cat_center[0] /= count2
        cat_center[1] /= count2
        cat_center[2] /= count2
        inactive_dist = dist_3d(act_mol_center, cat_center)
        active_dists.append(active_dist)
        inactive_dists.append(inactive_dist)
        differences.append(inactive_dist-active_dist)
    num_bins = 50
    #n, bins, patches = plt.hist(differences, num_bins, normed=1, facecolor='green', alpha=0.75)
    #n, bins, patches = plt.hist(differences, bins=[-30, 0, 30], normed=1, facecolor='green', alpha=0.75)
    #plt.xlabel('Differences between Active and Inactive Distances')
    #plt.ylabel('Probability')
    #plt.title('Histogram of Differences in Distance for Pos: ' + str(pos))
    ax = plt.figure().add_subplot(111)
    ax.set_xlim(0, 40)
    ax.set_ylim(-1, 1)
    y = 0.0
    # for d in differences:
        # plt.plot(d, y, 'ro', ms = 15, mfc='r')
    for a in active_dists:
        plt.plot(a, y+0.5, 'ro', ms = 15, alpha=0.1, mfc='r')
    for i in inactive_dists:
        plt.plot(i, y-0.5, 'bo', ms=15, alpha=0.1, mfc='b')
    plt.show()
    return differences, active_dists, inactive_dists

# Given a kinase for which PDBs for the active and inactive form exist within actives/complete/ and inactives/complete/, for each residue in the activation loop calculate the distance between that residue and the center of the catalytic loop for both the active and inactive states for all kinases listed in a file whose name is file_name. Compute p-values for a matched pairs test of the distances between the activation loop and the catalytic loop in the active and inactive states.
def findBestPos(file_name, endpoints_file, align_file, kinase_name):
    list_file = open(file_name, "r+")
    list_lines = list_file.readlines()
    kinases = []
    endpoint_info = extract_endpoints(endpoints_file)
    i = 0
    act_start = int(endpoint_info[kinase_name][len(endpoint_info[kinase_name])-2])
    act_end = int(endpoint_info[kinase_name][len(endpoint_info[kinase_name])-1])
    p_values = []
    for pos in range(act_start, act_end+1):
        while (i < len(list_lines)-1):
            line1 = list_lines[i]
            line2 = list_lines[i+1]
            line1_elements = line1.split(",")
            line2_elements = line2.split(",")
            if (line1_elements[1] == line2_elements[1]):
                to_add = []
                if (line1_elements[2] == "yes"):
                    to_add.append("actives/complete/" + line1_elements[1] + "_active.pdb")
                else:
                    to_add.append("actives/incomplete/" + line1_elements[1] + "_active.pdb")
                if (line2_elements[2] == "yes"):
                    to_add.append("inactives/complete/" + line2_elements[1] + "_inactive.pdb")
                else:
                    to_add.append("inactives/incomplete/" + line2_elements[1] + "_inactive.pdb")
                kinases.append(to_add)
                i += 2
            else:
                i += 1
        differences = []
        align_f = open(align_file, "r+")
        align_info = align_f.readlines()
        align_num, running_align_line = find_align_num(align_file, kinase_name, pos)
        active_dists = []
        inactive_dists = []
        for k in kinases:
            if (not os.path.isfile(k[0])) or (not os.path.isfile(k[1])): 
                continue
            active_file = open(k[0], "r+")
            inactive_file = open(k[1], "r+")
            active_lines = active_file.readlines()
            inactive_lines = inactive_file.readlines()
            k_name = ""
            if ("incomplete" in k[0]):
                k_name = k[0][ACTIVE_START_INCOMP:len(k[0])-ACTIVE_END]
            else:
                k_name = k[0][ACTIVE_START:len(k[0])-ACTIVE_END]
            kinase_endpoints = endpoint_info[k_name]
            cat_start = int(kinase_endpoints[len(kinase_endpoints)-4])
            cat_end = int(kinase_endpoints[len(kinase_endpoints)-3])
            i = 0
            active_dist = 0
            inactive_dist = 0
            equivalent_pos = -1
            act_mol_center = [0.0, 0.0, 0.0]
            cat_center = [0.0, 0.0, 0.0]
            count1 = 0
            count2 = 0
            for j in range(running_align_line, len(align_info)):
                split_line = align_info[j].split(" ")
                if (len(split_line) == 3):
                    if (split_line[0] != "align"):
                        if (split_line[0][1:] == k_name):
                            equivalent_pos = int(split_line[2][:len(split_line[2])-1])
                            break
            while (i < len(active_lines)):
                if ("ATOM" in active_lines[i]):
                    if (int(active_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) == equivalent_pos):
                        x_value = float(active_lines[i][POS_X_START:POS_X_END+1])
                        y_value = float(active_lines[i][POS_Y_START:POS_Y_END+1])
                        z_value = float(active_lines[i][POS_Z_START:POS_Z_END+1])
                        act_mol_center[0] += x_value
                        act_mol_center[1] += y_value
                        act_mol_center[2] += z_value
                        count1 += 1
                    elif (int(active_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= cat_start and int(active_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= cat_end):
                        x_value = float(active_lines[i][POS_X_START:POS_X_END+1])
                        y_value = float(active_lines[i][POS_Y_START:POS_Y_END+1])
                        z_value = float(active_lines[i][POS_Z_START:POS_Z_END+1])
                        cat_center[0] += x_value
                        cat_center[1] += y_value
                        cat_center[2] += z_value
                        count2 += 1
                i += 1
            act_mol_center[0] /= count1
            act_mol_center[1] /= count1
            act_mol_center[2] /= count1
            cat_center[0] /= count2
            cat_center[1] /= count2
            cat_center[2] /= count2
            active_dist = dist_3d(act_mol_center, cat_center)
            i = 0
            act_mol_center = [0.0, 0.0, 0.0]
            cat_center = [0.0, 0.0, 0.0]
            count1 = 0
            count2 = 0
            while (i < len(inactive_lines)):
                if ("ATOM" in inactive_lines[i]):
                    if (int(inactive_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1] ) == equivalent_pos):
                        x_value = float(inactive_lines[i][POS_X_START:POS_X_END+1])
                        y_value = float(inactive_lines[i][POS_Y_START:POS_Y_END+1])
                        z_value = float(inactive_lines[i][POS_Z_START:POS_Z_END+1])
                        act_mol_center[0] += x_value
                        act_mol_center[1] += y_value
                        act_mol_center[2] += z_value
                        count1 += 1
                    elif (int(inactive_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= cat_start and int(inactive_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= cat_end):
                        x_value = float(inactive_lines[i][POS_X_START:POS_X_END+1])
                        y_value = float(inactive_lines[i][POS_Y_START:POS_Y_END+1])
                        z_value = float(inactive_lines[i][POS_Z_START:POS_Z_END+1])
                        cat_center[0] += x_value
                        cat_center[1] += y_value
                        cat_center[2] += z_value
                        count2 += 1
                i += 1

            act_mol_center[0] /= count1
            act_mol_center[1] /= count1
            act_mol_center[2] /= count1
            cat_center[0] /= count2
            cat_center[1] /= count2
            cat_center[2] /= count2
            inactive_dist = dist_3d(act_mol_center, cat_center)
            active_dists.append(active_dist)
            inactive_dists.append(inactive_dist)
            differences.append(inactive_dist-active_dist)
        p_values.append([pos, spy.ttest_rel(active_dists, inactive_dists)[1]])
    p_values.sort(key=lambda x:x[1])
    return p_values

# Given the name of a file that contains a list of kinases for which both active and inactive PDBs are available, return a list of the names of all of these kinases.
def collect_key_kinases(file_name):
    lines = (open(file_name, "r+")).readlines()
    l = 1
    collection = []
    while l < len(lines)-1:
        kinase_name = lines[l].split(",")[1]
        next_kinase_name = lines[l+1].split(",")[1]
        inc = 1
        if (kinase_name == next_kinase_name):
            if ("active" in lines[l] and "inactive" in lines[l+1]):
                collection.append(kinase_name)
                inc = 2
        l += inc
    return collection

# Call the findBestPos() function for several different kinases (for which both active and inactive PDBs are available)
def find_good_align_nums(file_name, endpoints_file, align_file):
    key_kinases = collect_key_kinases(file_name)
    align_lines = []
    for k in key_kinases:
        print "Iteration " + str(k)
        print len(align_lines)
        print align_lines
        ordered_p_values = findBestPos(file_name, endpoints_file, align_file, k)
        t = 0
        while ordered_p_values[t][1] < 0.05:
            pair = find_align_num(align_file, k, ordered_p_values[t][0])
            if (pair[1] not in align_lines):
                align_lines.append(pair[1])
            t += 1
    return align_lines

# Compute a list whose size is final_size that contains those residues that are shared by the highest number of kinases at the same align numbers
def most_common_residues(align_file, final_size):
    align_lines = (open(align_file, "r+")).readlines()
    residue_pq = []
    heap_size = 0
    residue_dict = {}
    running_align_num = 0
    for l in range(0, len(align_lines)):
        if ("align" in align_lines[l]):
            for k in residue_dict.keys():
                if (heap_size < final_size):
                    heappush(residue_pq, (residue_dict[k], (k, running_align_num)))
                    heap_size += 1
                elif (residue_dict[k] > residue_pq[0][0]):
                    heappop(residue_pq)
                    heappush(residue_pq, (residue_dict[k], (k, running_align_num)))
            running_align_num = int(align_lines[l].split(" ")[2])
            residue_dict = {}
        else:
            line_str = align_lines[l].split(" ")
            residue_label = line_str[1]
            if (residue_label in residue_dict.keys()):
                residue_dict[residue_label] += 1
            else:
                residue_dict[residue_label] = 1
    output = []
    while (len(residue_pq) > 0):
        pq_element = heappop(residue_pq)
        output.append(pq_element)
    return output

# Given a residue, a kinase_name, an align_number, and a file with information about align_number, compute the number of the residue associated with that number within the given kinase or -1 if that kinase does not have the given residue at the specified align_number
def find_residue_num(align_file, align_num, residue, kinase_name):
    align_lines = (open(align_file, "r+")).readlines()
    running_align_num = 0
    for l in range(0, len(align_lines)):
        if ("align" in align_lines[l]):
            running_align_num = int(align_lines[l].split(" ")[2])
            if (running_align_num > align_num):
                break
        elif (running_align_num == align_num):
            if (align_lines[l].split(" ")[0][1:] == kinase_name and align_lines[l].split(" ")[1] == residue):
                return int(align_lines[l].split(" ")[2])
    return -1
