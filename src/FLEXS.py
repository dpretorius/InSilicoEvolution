# FLEXS.py
import logging
from utils import record_results
from colabfold_runner import run_colabfold_for_sequences
from structure_processing import compute_scores

import flexs
from flexs import baselines

# class DesignModel(flexs.Model):
#     def _fitness_function(self, sequences):
#         # Setup and run colabfold
#         output_dirs, design_number, final_output_path = run_colabfold_for_sequences(
#             sequences,
#             args.parent_dir,
#             args.input_dir,
#             args.output_dir,
#             args.final_output_dir,
#             args.num_repeats
#         )
#         logging.info("Processes finished MODEL")
#         
#         # Compute scores
#         beta_solenoid, plddt_scores, plddt_avg, loss_term = compute_scores(output_dirs)
#         logging.info("beta_solenoid scores: %s", beta_solenoid)
#         logging.info("plddt scores: %s", plddt_avg)
#         
#         # Record results
#         record_results("full", sequences, beta_solenoid, plddt_avg, loss_term, design_number, final_output_path)
#         return loss_term
# 
#     def train(self, *args, **kwargs):
#         pass
# 
class DesignModel(flexs.Model):
    def __init__(self, config, **kwargs):
        super().__init__(**kwargs)
        self.config = config  # Save configuration for use in methods

    def _fitness_function(self, sequences):
        output_dirs, design_number, final_output_path = run_colabfold_for_sequences(
            sequences,
            self.config.parent_dir,       # <-- Use self.config instead of args.
            self.config.input_dir,
            self.config.output_dir,
            self.config.final_output_dir,
            self.config.num_repeats
        )
        logging.info("Processes finished MODEL")
        beta_solenoid, plddt_scores, plddt_avg, loss_term = compute_scores(output_dirs, self.config.solenoid_type)
        logging.info("beta_solenoid scores: %s", beta_solenoid)
        logging.info("plddt scores: %s", plddt_avg)
        record_results("full", sequences, beta_solenoid, plddt_avg, loss_term, design_number, final_output_path)
        return loss_term

# class DesignLandscape(flexs.Landscape):
#     def _fitness_function(self, sequences):
#         # Setup and run colabfold
#         output_dirs, design_number, final_output_path = run_colabfold_for_sequences(
#             sequences,
#             args.parent_dir,
#             args.input_dir,
#             args.output_dir,
#             args.final_output_dir,
#             args.num_repeats
#         )
#         logging.info("Processes finished LANDSCAPE")
#         
#         # Compute scores
#         beta_solenoid, plddt_scores, plddt_avg, loss_term = compute_scores(output_dirs)
#         logging.info("beta_solenoid scores: %s", beta_solenoid)
#         logging.info("plddt scores: %s", plddt_avg)
#         
#         # Record results
#         record_results("best", sequences, beta_solenoid, plddt_avg, loss_term, design_number, final_output_path)
#         
#         # Additional post-processing if design criteria are met
#         if (beta_solenoid and beta_solenoid[0] > args.min_solenoid and 
#             plddt_scores and all(score > args.min_plddt for score in plddt_scores[0])):
#             results = [[sequences[i], beta_solenoid[i], plddt_avg[i], loss_term[i]] for i in range(len(sequences))]
#             handle_complete_design(output_dirs, design_number, final_output_path, args.output_dir, results)
#         return loss_term
        
class DesignLandscape(flexs.Landscape):
    def __init__(self, config, **kwargs):
        super().__init__(**kwargs)
        self.config = config

    def _fitness_function(self, sequences):
        output_dirs, design_number, final_output_path = run_colabfold_for_sequences(
            sequences,
            self.config.parent_dir,
            self.config.input_dir,
            self.config.output_dir,
            self.config.final_output_dir,
            self.config.num_repeats
        )
        logging.info("Processes finished LANDSCAPE")
        beta_solenoid, plddt_scores, plddt_avg, loss_term = compute_scores(output_dirs, self.config.solenoid_type)
        logging.info("beta_solenoid scores: %s", beta_solenoid)
        logging.info("plddt scores: %s", plddt_avg)
        record_results("best", sequences, beta_solenoid, plddt_avg, loss_term, design_number, final_output_path)
        if (beta_solenoid and beta_solenoid[0] > self.config.min_solenoid and 
            plddt_scores and all(score > self.config.min_plddt for score in plddt_scores[0])):
            results = [[sequences[i], beta_solenoid[i], plddt_avg[i], loss_term[i]] for i in range(len(sequences))]
            handle_complete_design(output_dirs, design_number, final_output_path, self.config.output_dir, results)
        return loss_term