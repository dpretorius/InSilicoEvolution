# InSilicoEvolution
Hallucination based protein design method with Alphafold2 as an oracle and SOLeNNoID discriminator network to produce solenoid proteins.

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
