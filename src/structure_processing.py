# structure_processing.py
import os
import json
import statistics
import numpy as np
from solenoid_detector_onehot import get_final_scores
from utils import compute_confidence_score

def find_average_plddt(outpur_dirs):
    normalised_plddt_list = []
    average_plddt_score = []
    for i in outpur_dirs:
        plddt_list = []
        for j in os.listdir(i):
            if os.path.isfile(os.path.join(i,j)) and 'scores.json' in j and 'unrelaxed_rank' in j:
                colabfold_output_file = os.path.join(i, j)
                input_file = open(colabfold_output_file, 'r')
                json_decode = json.load(input_file)

                plddt = json_decode['plddt']

                plddt_mu = np.mean(plddt)
                
                normalised_plddt = plddt_mu/100 # 0-1, 0 being the worst 1 being the best
                
                plddt_list.append(normalised_plddt)
                
        normalised_plddt_list.append(plddt_list)
        norm_plddt = statistics.mean(plddt_list)
        average_plddt_score.append(norm_plddt)
        
    return normalised_plddt_list, average_plddt_score                                     

# use solenoid detector 
def find_solenoid_probability(output_dirs, solenoid_type="beta"):
    combined_scores = []
    solenoid_map = {"beta": 1, "alphabeta": 2, "alpha": 3}
    desired_index = solenoid_map.get(solenoid_type, 1)
    
    for directory in output_dirs:
        individual_scores = []
        for filename in os.listdir(directory):
            if os.path.isfile(os.path.join(directory, filename)) and '.pdb' in filename and 'model' in filename:
                structure_path = os.path.join(directory, filename)
                final_scores = get_final_scores(structure_path)
                a = np.array(final_scores)
                one_hot = np.zeros_like(a)
                one_hot[np.arange(len(a)), a.argmax(axis=1)] = 1
                mean_scores = np.mean(one_hot, axis=0)[0:4]
                # Compute score: desired category score minus sum of the others
                score = mean_scores[desired_index] - sum(mean_scores[i] for i in range(4) if i != desired_index)
                individual_scores.append(score)
        combined_score = statistics.mean(individual_scores) if individual_scores else 0
        combined_scores.append(combined_score)
    return combined_scores

def compute_scores(output_dirs, solenoid_type="beta"):
    beta_solenoid_confidence = find_solenoid_probability(output_dirs, solenoid_type)
    plddt_scores, plddt_score = find_average_plddt(output_dirs)
    loss_term = compute_confidence_score(beta_solenoid_confidence, plddt_score)
    return beta_solenoid_confidence, plddt_scores, plddt_score, loss_term