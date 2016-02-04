from ML_Models import *
import random

surface_areas = {"A": 67, "R": 196, "N": 113, "D": 106, "C": 104, "Q": 144, "E": 138, "G": 0, "H": 151, "I": 140, "L": 137, "K": 167, "M": 160, "F": 175, "P": 105, "S": 80, "T": 102, "W": 217, "Y": 187, "V": 117}
volumes = {"A": 67, "R": 148, "N": 96, "D": 91, "C": 86, "Q": 114, "E": 109, "G": 48, "H": 118, "I": 124, "L": 124, "K": 135, "M": 124, "F": 135, "P": 90, "S": 73, "T": 93, "W": 163, "Y": 141, "V": 105}
polarities = {"A": 8.1, "K": 11.3, "D": 13.0, "P": 8.0, "E": 12.3, "W": 5.4, "I": 5.2, "L": 4.1, "N": 11.6, "F": 5.2, "Q": 10.5, "T": 8.6, "H": 10.4, "V": 5.9, "R": 10.5, "M": 5.7, "C": 5.5, "S": 9.2, "G": 9.0, "Y": 6.2}
significant_residues = [70518, 74845, 74572, 70059, 73367, 74426, 71862, 74278, 72249, 74713, 73955, 76274, 75279, 75844, 74142]

def compute_sphericity(amino):
    return float(np.pi*((6*float(volumes[amino]))**(2.0/3))/float(surface_areas[amino]))

def partition_dataset(filename):
    datafile = open(filename, "r+")
    datalines = datafile.readlines()
    training = []
    testing = []
    count = 0
    while (count < len(datalines)):
        if (len(training) >= len(datalines)/2):
            testing.append(count)
        elif (len(testing) >= len(datalines)/2):
            training.append(count)
        elif (random.random() >= 0.5):
            training.append(count)
        else:
            testing.append(count)
        count += 1
        return [training, testing] 
def compute_correlations(filename, align_file_name):
    datafile = open(filename, "r+")
    datalines = datafile.readlines()
    features = []
    pol = []
    pol_diffs = []
    pol_ratios = []
    spher = []
    spher_diffs = []
    spher_ratios = []
    residue_flags = []
    significant = []
    for r in range(0, len(significant_residues)):
        residue_flags.append([])
    classifications = []
    for i in range(0, len(datalines)):
        line = datalines[i].split(" ")
        kinase = line[0]
        print kinase + " " + str(i)
        wild_type = line[1]
        print wild_type
        res_num = int(line[2])
        mut_type = line[3]
        category = int(line[4])
        classifications.append(category)
        align_file_line_num = find_align_num(align_file_name, kinase, res_num)[1]
        flag = significant_residues.index(align_file_line_num) if align_file_line_num in significant_residues else -1
        if (flag != -1):
            residue_flags[flag].append(1)
            significant.append(1)
        else:
            significant.append(-1)
        for j in range(0, len(residue_flags)):
            if (j != flag):
                residue_flags[j].append(-1)
        if surface_areas[mut_type] != 0:
            spher.append(compute_sphericity(mut_type))
            if surface_areas[wild_type] != 0:
                spher_diffs.append(compute_sphericity(mut_type)-compute_sphericity(wild_type))
                spher_ratios.append(compute_sphericity(mut_type)/compute_sphericity(wild_type))
            else:
                spher_diffs.append(-1)
                spher_ratios.append(-1)
        else:
            spher.append(-1)
            spher_diffs.append(-1)
            spher_ratios.append(-1)
        pol.append(polarities[mut_type])
        pol_diffs.append(polarities[mut_type]-polarities[wild_type])
        pol_ratios.append(polarities[mut_type]/polarities[wild_type])
    features.append(classifications)
    features.append(pol)
    features.append(pol_diffs)
    features.append(pol_ratios)
    features.append(spher)
    features.append(spher_diffs)
    features.append(spher_ratios)
    for f in range(0, len(residue_flags)):
        features.append(residue_flags[f])
    features.append(significant)
    return np.corrcoef(features)
