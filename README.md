# Mutantsim

# MutantSim — Expected Protein Outcomes under Biased Polymerase Errors

> Given an mRNA and a biased base-substitution model per amplification round,
> estimate the protein-level consequences (unchanged AA, #nonsilent codons, premature/late STOP, moved start), 
> plus region-specific risks and synonymous recoding suggestions.

## What does this project do?
- Computes **expected** protein outcomes analytically using a 4×4 mutation matrix and `R` rounds.
- Optionally **simulates** many sequences to validate or explore edge cases.
- Reports region-specific risks and suggests synonymous codons to reduce (or increase) selected risks.

## Inputs and Outputs
**Input**
- mRNA sequence (FASTA or raw string; T is accepted and converted to U).
- 4×4 biased mutation matrix CSV (rows A,C,G,U → cols A,C,G,U) per round.
- Rounds `R` (integer).  
- Optional: region(s) of interest, downstream FASTA for homology scan.

**Output**
- Sequence length; original protein translation.
- Expected % unchanged protein; distribution over {k nonsilent codons}.
- % premature STOP; % delayed STOP; % moved start.
- ROI: % with ≥1 nonsilent; distribution over nonsilent count.
- (Optional) homology flags and synonymous recoding suggestions.
- HTML/Markdown report + CSV tables.

## Quickstart
```bash
git clone https://github.com/<you>/mutantsim.git
cd mutantsim
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest -q
mutantsim analyze --help
