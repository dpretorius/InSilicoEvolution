# InSilicoEvolution
Hallucination based protein design method with Alphafold2 as an oracle and SOLeNNoID discriminator network to produce solenoid proteins.

🚧 WIP 🚧

## Installation

```bash
git clone https://github.com/yourusername/InSilicoEvolution.git
cd InSilicoEvolution
```

## Running the program

```bash
python3 src/main.py \
  --parent_dir /path/to/output \
  --population_size 20 \
  --rounds 50 \
  --solenoid_type alphabeta
```

## Command-Line Arguments

| Argument                  | Description                                                  | Default                          |
|--------------------------|--------------------------------------------------------------|----------------------------------|
| `--parent_dir`           | Root directory for input/output                              | `.`                              |
| `--input_dir`            | Input FASTA folder name                                      | `in_silico_evolution_input`      |
| `--output_dir`           | ColabFold output folder name                                 | `in_silico_evolution_output`     |
| `--final_output_dir`     | Where final results are stored                               | `output_statistics`              |
| `--num_repeats`          | Repeats of the sequence in FASTA                             | `6`                              |
| `--population_size`      | Genetic algorithm population size                            | `10`                             |
| `--parent_strategy`      | Parent selection strategy                                    | `wright-fisher`                  |
| `--beta`                 | Mutation strength parameter                                  | `0.1`                            |
| `--children_proportion`  | Proportion of children per generation                        | `0.8`                            |
| `--rounds`               | Number of design rounds                                      | `30`                             |
| `--sequences_batch_size` | Number of sequences processed per batch                      | `1`                              |
| `--model_queries_per_batch` | Number of queries per generation                         | `30`                             |
| `--starting_sequence`    | Provide a starting sequence                                  | `""` (auto-generated)            |
| `--sequence_length`      | Length of generated starting sequence                        | `30`                             |
| `--min_solenoid`         | Threshold for solenoid confidence                            | `0.6`                            |
| `--min_plddt`            | Threshold for pLDDT confidence                               | `0.7`                            |
| `--solenoid_type`        | Solenoid class to target (`beta`, `alphabeta`, `alpha`)      | `beta`                           |

