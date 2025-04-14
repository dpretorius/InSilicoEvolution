# colabfold_runner.py
import os
from utils import round_counter_external, round_counter_internal

def execute_colabfold(input_dir, output_dir):
    # Run colabfold for a single input and output
    os.system(f'colabfold_batch {input_dir} {output_dir} --overwrite-existing-result --msa-mode single_sequence')

# write sequences into a fasta file and cmd args for colabfold
def prepare_fasta_files(seqs, input_path, output_path, num_repeats):
    num_rounds = round_counter_external()
    internal_rounds = round_counter_internal()
    design_number = num_rounds - internal_rounds
    count = 0
    output_dirs = []
    # For each sequence, create a FASTA file
    for index, sequence in enumerate(seqs, start=1):
        count += 1
        fasta_name = f'Round_{design_number}_{internal_rounds}_Sequence_{count}.fasta'
        name = f'Round_{design_number}_{internal_rounds}/Sequence_{count}'
        output_name = os.path.join(output_path, name)
        fasta_file = os.path.join(input_path, fasta_name)
        fasta_seq = sequence * num_repeats
        lines = [f'>Sequence_{count}', fasta_seq]
        with open(fasta_file, "w") as file1:
            file1.write('\n'.join(lines))
        output_dirs.append(output_name)
    return fasta_file, output_name, output_dirs, design_number
    
def run_colabfold_for_sequences(sequences, parent_dir, input_dir, output_dir, final_output_dir, num_repeats):
    input_path, out_path, final_out_path = setup_directories(parent_dir, input_dir, output_dir, final_output_dir)
    fasta_file, output_name, output_dirs, design_number = prepare_fasta_files(sequences, input_path, out_path, num_repeats)
    execute_colabfold(fasta_file, output_name)
    return output_dirs, design_number, final_out_path