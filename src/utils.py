# utils.py
import os
import csv
import statistics
from random import choice
import flexs.utils.sequence_utils as s_utils

# Alphabet
alphabet = s_utils.AAS

# Counters 
def round_counter_internal():
    round_counter_internal.counter += 1
    return round_counter_internal.counter
round_counter_internal.counter = 0

def round_counter_external():
    round_counter_external.counter += 1
    return round_counter_external.counter
round_counter_external.counter = 0

def setup_directories(parent_dir, input_dir, output_dir, final_output_dir):
    input_path = os.path.join(parent_dir, input_dir)
    output_path = os.path.join(parent_dir, output_dir)
    final_output_path = os.path.join(output_path, final_output_dir)
    for path in (input_path, output_path, final_output_path):
        os.makedirs(path, exist_ok=True)
    return input_path, output_path, final_output_path
    
#random sequence generator
def generate_initial_sequence(length_of_seq):
    protein_seq=""
    for count in range(length_of_seq):
        protein_seq+=choice(alphabet)
    return protein_seq
    
def append_results_to_csv(filename, results, design_number, output_dir):
    for i in range(0,len(results)):
        results[i].insert(0, design_number)
    with open(f"{output_dir}/{filename}.csv", 'a+') as f:
        writer = csv.writer(f)
        for i in range (len(results)):
            writer.writerow(results[i])
    return 
    
# combine difference metrics for loss function 
def compute_confidence_score(combined_scores, tm_scores):
    confidence_scores = []
    for (item1, item2) in zip(combined_scores, tm_scores):
        confidence_scores.append(item1+item2)
    return confidence_scores
    
def record_results(filename, sequences, beta_solenoid, plddt, loss, design_number, output_dir):
    results = [[sequences[i], beta_solenoid[i], plddt[i], loss[i]] for i in range(len(sequences))]
    append_results_to_csv(filename, results, design_number, output_dir)
    
def handle_complete_design(output_dirs, design_number, final_output_path, output_path, results):
    complete_folder = os.path.join(final_output_path, f'Complete_Design_Round_{design_number}')
    os.makedirs(complete_folder, exist_ok=True)
    for file_name in os.listdir(output_dirs[0]):
        source = os.path.join(output_dirs[0], file_name)
        destination = os.path.join(complete_folder, file_name)
        if os.path.isfile(source):
            shutil.copy(source, destination)
    # Record the top result
    append_results_to_csv("complete", [results[0]], design_number, final_output_path)
    
    # Zip files and clean up
    zip_name = os.path.join(output_path, f'Design_Round_{design_number}')
    os.makedirs(zip_name, exist_ok=True)
    for folder in os.listdir(output_path):
        if "Round" in folder:
            src = os.path.join(output_path, folder)
            shutil.move(src, zip_name)
    shutil.make_archive(zip_name, 'zip', os.path.join(output_path, os.path.basename(zip_name)))
    shutil.rmtree(zip_name)
    sys.exit()