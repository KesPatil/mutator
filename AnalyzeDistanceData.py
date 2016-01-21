import math
import scipy
import scipy.constants
import scipy.stats as spy
import numpy as np
import matplotlib.pyplot as plt
import random
import os

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
AMINO_START = 18-1
AMINO_END = 20-1
ATOM_START = 13-1
ATOM_END = 16-1


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

def get_training_testing_sets(ac_fn, inac_fn, flag):
    ac_file = open(ac_fn, "r+")
    inac_file = open(inac_fn, "r+")
    ac_lines = ac_file.readlines()
    inac_lines = inac_file.readlines()
    training = []
    testing = []
    for acl in ac_lines:
        if (flag == 1):
            training.append([1, acl[:len(acl)-1]])
        elif (len(training) == len(ac_lines)/2):
            testing.append([1, acl[:len(acl)-1]])
        elif (len(testing) == len(ac_lines)/2):
            training.append([1, acl[:len(acl)-1]])
        elif (random.random() < 0.5):
            training.append([1, acl[:len(acl)-1]])
        else:
            testing.append([1, acl[:len(acl)-1]])
    for inacl in inac_lines:
        if (flag == 1):
            training.append([0, inacl[:len(inacl)-1]])
        elif (len(training) - len(ac_lines)/2 >= len(inac_lines)/2):
            print "no"
            testing.append([0, inacl[:len(inacl)-1]])
        elif (len(testing) - len(ac_lines)/2 >= len(inac_lines)/2):
            training.append([0, inacl[:len(inacl)-1]])
        elif (random.random() < 0.5):
            training.append([0, inacl[:len(inacl)-1]])
        else:
            print "no"
            testing.append([0, inacl[:len(inacl)-1]])
    return [training, testing]

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

def dist_3d(xyz1, xyz2):
    return math.sqrt((xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 + (xyz1[2]-xyz2[2])**2)

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

def find_next_amino(pdb_lines, currIndex):
    if ("ATOM" not in pdb_lines[currIndex]):
        return ["", -1]
    currAmino = pdb_lines[currIndex][AMINO_START:AMINO_END+1]
    while (pdb_lines[currIndex][AMINO_START:AMINO_END+1] == currAmino):
        currIndex += 1
    return [pdb_lines[currIndex][AMINO_START:AMINO_END+1], currIndex]

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

def electric_field(pos1, charge2, pos2):
    r_vec = pos2 - pos1
    mult_constant = 1.0/(4*scipy.constants.pi*scipy.constants.epsilon_0)
    return mult_constant*charge2*r_vec/(np.linalg.norm(r_vec)**3)

def electric_potential(pos1, charge2, pos2):
    r_vec = pos2-pos1
    mult_constant = 1.0/(4*scipy.constants.pi*scipy.constants.epsilon_0)
    return mult_constant*charge2/np.linalg.norm(r_vec)

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
    print p_values

