from AnalyzeDistanceData import *
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.svm import SVC
from Bio.PDB import *
from Bio.PDB.HSExposure import *
from Bio.PDB.DSSP import *

PATH_TO_DSSP = "/usr/local/Cellar/dssp/2.2.1/bin/mkdssp"
PHI_HELIX = -64
PSI_HELIX = -41
PHI_PSI_DEVIATION = 7
amino_code_dict = {"F": "PHE", "L": "LEU", "S": "SER", "Y": "TYR", "X": "TER", "C": "CYS", "W": "TRP", "P": "PRO", "H": "HIS", "Q": "GLN", "R": "ARG", "I": "ILE", "M": "MET", "T": "THR", "N": "ASN", "K": "LYS", "V": "VAL", "A": "ALA", "D": "ASP", "E": "GLU", "G": "GLY"}

# Old function for constructing datasets
def get_data_sets_old(data_files, endpoints_file, align_file_name):
    dataset = []
    classifications = []
    align_file = open(align_file_name, "r+")
    align_lines = align_file.readlines()
    # ABL1 388?
    # align_num = 615
    # ABL1 384
    #align_num = 611
    # ABL1 385
    # align_num = 612
    #pos_list = [385, 397, 395, 384, 392]
    #pos_kinase_names = ["ABL1", "ABL1", "ABL1", "ABL1", "ABL1"]
    pos_list = [385, 397, 395]
    pos_kinase_names = ["ABL1", "ABL1", "ABL1"]
    pos_align_lines = [70518, 74845, 74572, 70059, 73367, 74426, 71862, 74278, 72249, 74713, 73955, 76274, 75279, 75844, 74142]
    align_num_list = []
    #for p in range(0, len(pos_list)):
        #align_num_list.append(find_align_num(align_file_name, pos_kinase_names[p], pos_list[p])[0])
    for p in range(0, len(pos_align_lines)):
        align_num_list.append(int(align_lines[p].split(" ")[2][:len(align_lines[p].split(" ")[2])-1]))
    endpoint_info = extract_endpoints(endpoints_file)
    for t in data_files:
        data = []
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
        alphac_start = int(kinase_endpoints[len(kinase_endpoints)-6])
        alphac_end = int(kinase_endpoints[len(kinase_endpoints)-5])
        first_residue = int(kinase_endpoints[1])

        parser = PDBParser()
        structure = parser.get_structure(kinase_name, t[1])
        model = structure[0]
        dssp = DSSP(model, t[1], dssp=PATH_TO_DSSP)
        dssp_list = list(dssp)
        phi_deviations = []
        psi_deviations = []
        for ind in range(alphac_start, alphac_end+1):
            dssp_list_index = ind-first_residue+1
            phi = dssp_list[dssp_list_index][4]
            psi = dssp_list[dssp_list_index][5]
            phi_deviations.append(np.sqrt((phi-PHI_HELIX)**2))
            psi_deviations.append(np.sqrt((psi-PSI_HELIX)**2))
        data.append(np.mean(phi_deviations))
        data.append(np.mean(psi_deviations))

        for c in range(0, len(align_num_list)):
            align_num = align_num_list[c]
            running_align_num = -1
            equivalent_pos = -1
            for j in range(0, len(align_lines)):
                split_line = align_lines[j].split(" ")
                if (len(split_line) == 3):
                    if (split_line[0] == "align"):
                        running_align_num = int(split_line[2][:len(split_line[2])-1])
                    else:
                        if (split_line[0][1:] == kinase_name and running_align_num == align_num):
                            equivalent_pos = int(split_line[2][:len(split_line[2])-1])
                            break
            cat_center_pos = [0.0, 0.0, 0.0]
            act_mol_center = [0.0, 0.0, 0.0]
            count_act_mol = 0
            i = 0
            while (i < len(pdb_lines)):
                if ("ATOM" in pdb_lines[i]):
                    if ((int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= cat_start) and (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= cat_end)):
                        x_value = float(pdb_lines[i][POS_X_START:POS_X_END+1])
                        y_value = float(pdb_lines[i][POS_Y_START:POS_Y_END+1])
                        z_value = float(pdb_lines[i][POS_Z_START:POS_Z_END+1])
                        cat_center_pos = [cat_center_pos[0]+x_value, cat_center_pos[1]+y_value, cat_center_pos[2]+z_value]
                    elif (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) == equivalent_pos):
                        x_value = float(pdb_lines[i][POS_X_START:POS_X_END+1])
                        y_value = float(pdb_lines[i][POS_Y_START:POS_Y_END+1])
                        z_value = float(pdb_lines[i][POS_Z_START:POS_Z_END+1])
                        act_mol_center[0] += x_value
                        act_mol_center[1] += y_value
                        act_mol_center[2] += z_value
                        count_act_mol += 1
                i += 1
            cat_center_pos = np.array(cat_center_pos)/(cat_end-cat_start+1)
            if (count_act_mol > 0):
                act_mol_center[0] /= count_act_mol
                act_mol_center[1] /= count_act_mol
                act_mol_center[2] /= count_act_mol
                #data.append(dist_3d(cat_center_pos, act_mol_center))
            #else:
                #data.append(0)
        i = 0
        flag = 0
        flag2 = 0
        cat_field = np.array([0.0, 0.0, 0.0])
        min_distances = []
        X = []
        Z = []
        X_alpha = []
        Z_alpha = []
        while (i < len(pdb_lines)):
            if ("ATOM" in pdb_lines[i]):
                if (flag2 == 0):
                    flag2 = 1
                    ind = i
                if ((int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= act_start) and (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= act_end)):
                    flag = 1
                    if (pdb_lines[i][ATOM_START:ATOM_END+1].strip(' ') == "CA"):
                        the_dist_list = dist_list(pdb_lines, i, ind, cat_start, cat_end)
                        min_distances.append(min(the_dist_list))
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
                    X.append([x_value, y_value, 1])
                    Z.append(z_value)
                elif (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= alphac_start) and (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= alphac_end):
                    x_value = float(pdb_lines[i][POS_X_START:POS_X_END+1])
                    y_value = float(pdb_lines[i][POS_Y_START:POS_Y_END+1])
                    z_value = float(pdb_lines[i][POS_Z_START:POS_Z_END+1])
                    X_alpha.append([x_value, y_value, 1])
                    Z_alpha.append(z_value)
            i += 1
        data.append(sum(min_distances)/len(min_distances))
        #dataset.append([np.linalg.norm(cat_field)])
        #dataset[len(dataset)-1].append(np.linalg.norm(cat_field))
        X = np.array(X)
        Z = np.array(Z)
        least_squares = np.linalg.lstsq(X, Z)
        #dataset[len(dataset)-1].append(np.mean(least_squares[1]))
        least_squares_alpha = np.linalg.lstsq(X_alpha, Z_alpha)
        #data.append(np.mean(least_squares_alpha[1]))
        dataset.append(np.array(data))
        classifications.append(t[0])
    return [dataset, classifications]

# Construct either a dataset and a list of associated classifications (active or inactive) with various features to distinguish between active and inactive kinases.
def get_data_sets(data_files, endpoints_file, align_file_name):
    dataset = []
    classifications = []
    align_file = open(align_file_name, "r+")
    align_lines = align_file.readlines()
    pos_align_lines = [70518, 74845, 74572, 70059, 73367, 74426, 71862, 74278, 72249, 74713, 73955, 76274, 75279, 75844, 74142]
    common_residues = most_common_residues(align_file_name, 10)
    align_num_list = []
    for p in range(0, len(pos_align_lines)):
        align_num_list.append(int(align_lines[p].split(" ")[2][:len(align_lines[p].split(" ")[2])-1]))
    endpoint_info = extract_endpoints(endpoints_file)
    actCount = 0
    inactCount = 0
    for t_count in data_files:
        if (t_count[0] == 1 and "incomplete" in t_count[1]):
            kinase_name = t_count[1][ACTIVE_START_INCOMP:t_count[1].index("_")]
        elif (t_count[0] == 1):
            kinase_name = t_count[1][ACTIVE_START:t_count[1].index("_")]
        elif (t_count[0] == 0 and "incomplete" in t_count[1]):
            kinase_name = t_count[1][INACTIVE_START_INCOMP:t_count[1].index("_")]
        else:
            kinase_name = t_count[1][INACTIVE_START:t_count[1].index("_")]
        if (t_count[0] == 1 and kinase_name in endpoint_info):
            actCount += 1
        elif(t_count[0] == 0 and kinase_name in endpoint_info):
            inactCount += 1
    print actCount
    print inactCount
    actCount = min(actCount, inactCount)
    inactCount = actCount
    for t in data_files:
        data = []
        print t[1]
        pdb_file = open(t[1], "r+")
        pdb_lines = pdb_file.readlines()
        kinase_name = ""
        if (t[0] == 1 and "incomplete" in t[1]):
            kinase_name = t[1][ACTIVE_START_INCOMP:t[1].index("_")]
        elif (t[0] == 1):
            kinase_name = t[1][ACTIVE_START:t[1].index("_")]
        elif (t[0] == 0 and "incomplete" in t[1]):
            kinase_name = t[1][INACTIVE_START_INCOMP:t[1].index("_")]
        else:
            kinase_name = t[1][INACTIVE_START:t[1].index("_")]
        if (kinase_name not in endpoint_info):
            print kinase_name + " not in endpoint_info"
            continue
        kinase_endpoints = endpoint_info[kinase_name]
        cat_start = int(kinase_endpoints[len(kinase_endpoints)-4])
        cat_end = int(kinase_endpoints[len(kinase_endpoints)-3])
        act_start = int(kinase_endpoints[len(kinase_endpoints)-2])
        act_end = int(kinase_endpoints[len(kinase_endpoints)-1])
        alphac_start = int(kinase_endpoints[len(kinase_endpoints)-6])
        alphac_end = int(kinase_endpoints[len(kinase_endpoints)-5])
        ploop_start = int(kinase_endpoints[len(kinase_endpoints)-8])
        ploop_end = int(kinase_endpoints[len(kinase_endpoints)-7])
        first_residue = int(kinase_endpoints[0])

        parser = PDBParser()
        structure = parser.get_structure(kinase_name, t[1])
        model = structure[0]
        dssp = DSSP(model, t[1], dssp=PATH_TO_DSSP)
        dssp_list = list(dssp)
        phi_deviations = []
        psi_deviations = []
        act_accsurfarea = []
        cat_accsurfarea = []
        for ind in range(alphac_start, alphac_end+1):
            dssp_list_index = ind-first_residue+1
            phi = dssp_list[dssp_list_index][4]
            psi = dssp_list[dssp_list_index][5]
            phi_deviations.append((phi-PHI_HELIX))
            psi_deviations.append((psi-PSI_HELIX))
        data.append(np.mean(phi_deviations))
        data.append(np.mean(psi_deviations))
        for ind in range(act_start, act_end+1):
            dssp_list_index = ind-first_residue+1
            rel_acc = dssp_list[dssp_list_index][3]
            res_code = dssp_list[dssp_list_index][1]
            act_accsurfarea.append(rel_acc*MAX_ACC[amino_code_dict[res_code]])
        data.append(np.mean(act_accsurfarea))
#        for ind in range(cat_start, cat_end+1):
#            dssp_list_index = ind-first_residue+1
#            rel_acc = dssp_list[dssp_list_index][3]
#            res_code = dssp_list[dssp_list_index][1]
#            cat_accsurfarea.append(rel_acc*MAX_ACC[amino_code_dict[res_code]])
#        data.append(np.mean(cat_accsurfarea))

        residues = Selection.unfold_entities(structure, 'R')
        for res in common_residues:
            res_num = find_residue_num(align_file_name, res[1][1], res[1][0], kinase_name)-first_residue
            for res2 in common_residues:
                res2_num = find_residue_num(align_file_name, res2[1][1], res2[1][0], kinase_name)-first_residue
                if (res_num > -1 and residues[res_num].has_id('CA') and res2_num > -1 and residues[res2_num].has_id('CA') and (res[1][1] != res2[1][1] or res[1][0] != res2[1][0])):
                    data.append(residues[res_num]['CA']-residues[res2_num]['CA'])
                    origin = Vector((0, 0, 0))
                    data.append(calc_angle(origin, residues[res_num]['CA'].get_vector(), residues[res2_num]['CA'].get_vector()))
                else:
                    data.append(0)
                    data.append(0)

        for c in range(0, len(align_num_list)):
            #break
            align_num = align_num_list[c]
            running_align_num = -1
            equivalent_pos = -1
            for j in range(0, len(align_lines)):
                split_line = align_lines[j].split(" ")
                if (len(split_line) == 3):
                    if (split_line[0] == "align"):
                        running_align_num = int(split_line[2][:len(split_line[2])-1])
                    else:
                        if (split_line[0][1:] == kinase_name and running_align_num == align_num):
                            equivalent_pos = int(split_line[2][:len(split_line[2])-1])
                            break
            cat_center_pos = [0.0, 0.0, 0.0]
            act_mol_center = [0.0, 0.0, 0.0]
            count_act_mol = 0
            i = 0
            while (i < len(pdb_lines)):
                if ("ATOM" in pdb_lines[i]):
                    if ((int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= cat_start) and (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= cat_end)):
                        x_value = float(pdb_lines[i][POS_X_START:POS_X_END+1])
                        y_value = float(pdb_lines[i][POS_Y_START:POS_Y_END+1])
                        z_value = float(pdb_lines[i][POS_Z_START:POS_Z_END+1])
                        cat_center_pos = [cat_center_pos[0]+x_value, cat_center_pos[1]+y_value, cat_center_pos[2]+z_value]
                    elif (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) == equivalent_pos):
                        x_value = float(pdb_lines[i][POS_X_START:POS_X_END+1])
                        y_value = float(pdb_lines[i][POS_Y_START:POS_Y_END+1])
                        z_value = float(pdb_lines[i][POS_Z_START:POS_Z_END+1])
                        act_mol_center[0] += x_value
                        act_mol_center[1] += y_value
                        act_mol_center[2] += z_value
                        count_act_mol += 1
                i += 1
            cat_center_pos = np.array(cat_center_pos)/(cat_end-cat_start+1)
            if (count_act_mol > 0):
                act_mol_center[0] /= count_act_mol
                act_mol_center[1] /= count_act_mol
                act_mol_center[2] /= count_act_mol
                data.append(dist_3d(cat_center_pos, act_mol_center))
            else:
                data.append(0)
        i = 0
        flag = 0
        flag2 = 0
        cat_field = np.array([0.0, 0.0, 0.0])
        min_distances = []
        X = []
        Z = []
        X_alpha = []
        Z_alpha = []
        while (i < len(pdb_lines)):
            if ("ATOM" in pdb_lines[i]):
                if (flag2 == 0):
                    flag2 = 1
                    ind = i
                if ((int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= act_start) and (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= act_end)):
                    flag = 1
                    if (pdb_lines[i][ATOM_START:ATOM_END+1].strip(' ') == "CA"):
                        x_value = float(pdb_lines[i][POS_X_START:POS_X_END+1])
                        y_value = float(pdb_lines[i][POS_Y_START:POS_Y_END+1])
                        z_value = float(pdb_lines[i][POS_Z_START:POS_Z_END+1])
                        act_pos = np.array([x_value, y_value, z_value])
                        X.append([x_value, y_value, 1])
                        Z.append(z_value)
                elif (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) >= alphac_start) and (int(pdb_lines[i][RESIDUE_NUM_START:RESIDUE_NUM_END+1]) <= alphac_end):
                    if (pdb_lines[i][ATOM_START:ATOM_END+1].strip(' ') == "CA"):
                        x_value = float(pdb_lines[i][POS_X_START:POS_X_END+1])
                        y_value = float(pdb_lines[i][POS_Y_START:POS_Y_END+1])
                        z_value = float(pdb_lines[i][POS_Z_START:POS_Z_END+1])
                        X_alpha.append([x_value, y_value, 1])
                        Z_alpha.append(z_value)
            i += 1
        act_cat_distances = two_region_distances(kinase_name, pdb_lines, [act_start, act_end], [cat_start, cat_end])
        min_act_cat_dist = min_array(act_cat_distances)
        data.append(sum(min_act_cat_dist)/len(min_act_cat_dist))
        alpha_cat_distances = two_region_distances(kinase_name, pdb_lines, [alphac_start, alphac_end], [cat_start, cat_end])
        min_alpha_cat_dist = min_array(alpha_cat_distances)
        data.append(sum(min_alpha_cat_dist)/len(min_alpha_cat_dist))
        max_alpha_cat_dist = max_array(alpha_cat_distances)
        data.append(sum(max_alpha_cat_dist)/len(max_alpha_cat_dist))
        mean_alpha_cat_dist = mean_array(alpha_cat_distances)
        data.append(sum(mean_alpha_cat_dist)/len(mean_alpha_cat_dist))
        alpha_ploop_distances = two_region_distances(kinase_name, pdb_lines, [alphac_start, alphac_end], [ploop_start, ploop_end])
        min_alpha_ploop_dist = min_array(alpha_ploop_distances)
        max_alpha_ploop_dist = max_array(alpha_ploop_distances)
        mean_alpha_ploop_dist = mean_array(alpha_ploop_distances)
        data.append(sum(min_alpha_ploop_dist)/len(min_alpha_ploop_dist))
        data.append(sum(max_alpha_ploop_dist)/len(max_alpha_ploop_dist))
        data.append(sum(mean_alpha_ploop_dist)/len(mean_alpha_ploop_dist))
        act_ploop_distances = two_region_distances(kinase_name, pdb_lines, [act_start, act_end], [ploop_start, ploop_end])
        min_act_ploop_dist = min_array(act_ploop_distances)
        max_act_ploop_dist = max_array(act_ploop_distances)
        mean_act_ploop_dist = mean_array(act_ploop_distances)
        data.append(sum(min_act_ploop_dist)/len(min_act_ploop_dist))
        data.append(sum(max_act_ploop_dist)/len(max_act_ploop_dist))
        data.append(sum(mean_act_ploop_dist)/len(mean_act_ploop_dist))
        cat_ploop_distances = two_region_distances(kinase_name, pdb_lines, [cat_start, cat_end], [ploop_start, ploop_end])
        min_cat_ploop_dist = min_array(cat_ploop_distances)
        max_cat_ploop_dist = max_array(cat_ploop_distances)
        mean_cat_ploop_dist = mean_array(cat_ploop_distances)
        data.append(sum(min_cat_ploop_dist)/len(min_cat_ploop_dist))
        data.append(sum(max_cat_ploop_dist)/len(max_cat_ploop_dist))
        data.append(sum(mean_cat_ploop_dist)/len(mean_cat_ploop_dist))
        act_alpha_distances = two_region_distances(kinase_name, pdb_lines, [act_start, act_end], [alphac_start, alphac_end])
        min_act_alpha_dist = min_array(act_alpha_distances)
        max_act_alpha_dist = max_array(act_alpha_distances)
        mean_act_alpha_dist = mean_array(act_alpha_distances)
        data.append(sum(min_act_alpha_dist)/len(min_act_alpha_dist))
        data.append(sum(max_act_alpha_dist)/len(max_act_alpha_dist))
        data.append(sum(mean_act_alpha_dist)/len(mean_act_alpha_dist))
        data.append(min(min_act_cat_dist))
        data.append(min(min_act_ploop_dist))
        data.append(min(min_alpha_cat_dist))
        data.append(min(min_alpha_ploop_dist))
        data.append(min(min_cat_ploop_dist))
        X = np.array(X)
        Z = np.array(Z)
        least_squares = np.linalg.lstsq(X, Z)
        data.append(np.mean(least_squares[1]))
        least_squares_alpha = np.linalg.lstsq(X_alpha, Z_alpha)
        data.append(np.sqrt(np.mean(least_squares_alpha[1]**2)))
        if(True):
        #if ((t[0] == 1 and actCount > 0) or (t[0] == 0 and inactCount > 0)):
            dataset.append(np.array(data))
            classifications.append(t[0])
            if (t[0] == 1):
                actCount -= 1
            elif (t[0] == 0):
                inactCount -= 1
    print actCount
    print inactCount
    return [dataset, classifications]

# Compute a logistic regression to classify kinases based on the features in get_data_sets()
def kinase_logreg(training, testing, endpoints_file, align_file_name):
    [dataset, classifications] = get_data_sets(training, endpoints_file, align_file_name)
    logreg_model = LogisticRegression()
    logreg_model.fit(dataset, classifications)
    training_predictions = logreg_model.predict(dataset)
    print "Training set results:"
    print metrics.classification_report(classifications, training_predictions)
    print metrics.confusion_matrix(classifications, training_predictions)
    print metrics.roc_auc_score(classifications, training_predictions)
    if (len(testing) == 0):
        return [logreg_model.coef_, logreg_model.intercept_]
    print "Testing set results:"
    [testset, test_classifications] = get_data_sets(testing, endpoints_file, align_file_name)
    testing_predictions = logreg_model.predict(testset)
    print metrics.classification_report(test_classifications, testing_predictions)
    print metrics.confusion_matrix(test_classifications, testing_predictions)
    print metrics.roc_auc_score(test_classifications, testing_predictions)
    return [logreg_model.coef_, logreg_model.intercept_]

# Construct an SVM model to classify kinases as active or inactive based on the features in get_data_sets()
def kinase_svm(training, testing, endpoints_file, align_file_name, C_val=1.0):
    [dataset, classifications] = get_data_sets(training, endpoints_file, align_file_name)
    [testset, test_classifications] = get_data_sets(testing, endpoints_file, align_file_name)
    #svm_model = SVC(class_weight='balanced')
    svm_model = SVC(kernel="linear", C=C_val)
    svm_model.fit(dataset, classifications)
    training_predictions = svm_model.predict(dataset)
    print "Training set results:"
    print metrics.classification_report(classifications, training_predictions)
    print metrics.roc_auc_score(classifications, training_predictions)
    if (len(testing) > 0):
        print "Testing set results:"
        testing_predictions = svm_model.predict(testset)
        print metrics.classification_report(test_classifications, testing_predictions)
        print metrics.roc_auc_score(test_classifications, testing_predictions)
        return metrics.roc_auc_score(test_classifications, testing_predictions)
    return metrics.roc_auc_score(classifications, training_predictions)
