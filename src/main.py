# main.py
import argparse
import logging
from utils import generate_initial_sequence, round_counter_external
from FLEXS import DesignModel, DesignLandscape
from flexs import baselines
import flexs.utils.sequence_utils as s_utils

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Protein design pipeline")
    parser.add_argument("--parent_dir", type=str, default=".",
                        help="Parent directory path")
    parser.add_argument("--input_dir", type=str, default="in_silico_evolution_input",
                        help="Input directory name")
    parser.add_argument("--output_dir", type=str, default="in_silico_evolution_output",
                        help="Output directory name")
    parser.add_argument("--final_output_dir", type=str, default="output_statistics",
                        help="Final output directory name")
    parser.add_argument("--num_repeats", type=int, default=6,
                        help="Number of repeats for each sequence")
    parser.add_argument("--population_size", type=int, default=10,
                        help="Population size for the genetic algorithm")
    parser.add_argument("--parent_strategy", type=str, default="wright-fisher",
                        help="Parent selection strategy")
    parser.add_argument("--beta", type=float, default=0.1,
                        help="Beta value for the genetic algorithm")
    parser.add_argument("--children_proportion", type=float, default=0.8,
                        help="Children proportion for the genetic algorithm")
    parser.add_argument("--rounds", type=int, default=30,
                        help="Number of rounds for the genetic algorithm")
    parser.add_argument("--sequences_batch_size", type=int, default=1,
                        help="Sequences batch size (max GPU capacity)")
    parser.add_argument("--model_queries_per_batch", type=int, default=30,
                        help="Model queries per batch")
    parser.add_argument("--starting_sequence", type=str, default="",
                        help="Starting sequence; if empty, one will be generated")
    parser.add_argument("--sequence_length", type=int, default=30,
                        help="Length of the starting sequence if not provided")
    parser.add_argument("--min_solenoid", type=float, default=0.6,
                        help="Minimum solenoid score required for success")
    parser.add_argument("--min_plddt", type=float, default=0.7,
                        help="Minimum pLDDT score required for success")
    parser.add_argument("--solenoid_type", type=str, default="beta",
                        choices=["beta", "alphabeta", "alpha"],
                        help="Type of solenoid to evaluate (beta, alphabeta, alpha).")
    
    args = parser.parse_args()
    
    if not args.starting_sequence:
        args.starting_sequence = generate_initial_sequence(args.sequence_length)
    
    config = args
    
    # Instantiate model objects with configuration passed in their constructor
    design_model = DesignModel(config=config, name="DesignModel")
    design_landscape = DesignLandscape(config=config, name="DesignLandscape")
    
    def test_adalead():
        explorer = baselines.explorers.GeneticAlgorithm(
            model=design_model,
            population_size=config.population_size,
            parent_selection_strategy=config.parent_strategy,
            beta=config.beta,
            children_proportion=config.children_proportion,
            rounds=config.rounds,
            sequences_batch_size=config.sequences_batch_size,
            model_queries_per_batch=config.model_queries_per_batch,
            starting_sequence=config.starting_sequence,
            alphabet=s_utils.AAS,
        )
        adalead_sequences, metadata = explorer.run(design_landscape)
        return adalead_sequences, metadata

    for i in range(1, config.rounds + 1):
        round_counter_external()
        try:
            adalead_sequences, metadata = test_adalead()
        except Exception as e:
            logging.error("Error in iteration %d: %s", i, e)
            continue
