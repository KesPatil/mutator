from AnalyzeDistanceData import *
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.svm import SVC

def get_data_sets(data_files, endpoints_file):
    # Dataset: min_distance, field magnitude, planar similarity
    dataset = []
    classifications = []
    endpoint_info = extract_endpoints(endpoints_file)
    for t in data_files:
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
        dataset.append([sum(min_distances)/len(min_distances)])
        #dataset.append([np.linalg.norm(cat_field)])
        #dataset[len(dataset)-1].append(np.linalg.norm(cat_field))
        X = np.array(X)
        Z = np.array(Z)
        least_squares = np.linalg.lstsq(X, Z)
        dataset[len(dataset)-1].append(np.mean(least_squares[1]))
        least_squares_alpha = np.linalg.lstsq(X_alpha, Z_alpha)
        dataset[len(dataset)-1].append(np.mean(least_squares_alpha[1]))
        classifications.append(t[0])
    return [dataset, classifications]


def kinase_logreg(training, testing, endpoints_file):
    [dataset, classifications] = get_data_sets(training, endpoints_file)
    logreg_model = LogisticRegression()
    logreg_model.fit(dataset, classifications)
    training_predictions = logreg_model.predict(dataset)
    print "Training set results:"
    print metrics.classification_report(classifications, training_predictions)
    print metrics.confusion_matrix(classifications, training_predictions)
    print "Testing set results:"
    [testset, test_classifications] = get_data_sets(testing, endpoints_file)
    testing_predictions = logreg_model.predict(testset)
    print metrics.classification_report(test_classifications, testing_predictions)
    print metrics.confusion_matrix(test_classifications, testing_predictions)

def kinase_svm(training, testing, endpoints_file):
    [dataset, classifications] = get_data_sets(training, endpoints_file)
    [testset, test_classifications] = get_data_sets(testing, endpoints_file)
    svm_model = SVC()
    svm_model.fit(dataset, classifications)
    training_predictions = svm_model.predict(dataset)
    testing_predictions = svm_model.predict(testset)
    print "Training set results:"
    print metrics.classification_report(classifications, training_predictions)
    print "Testing set results:"
    print metrics.classification_report(test_classifications, testing_predictions)
