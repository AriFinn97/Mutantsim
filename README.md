# 🧬 MutantSim — Protein Outcomes Under Biased Polymerase Errors

## 📖 What does this project do?
MutantSim is a Python tool that takes an **mRNA coding sequence (CDS)** and a **mutation bias matrix** for an error-prone polymerase, and estimates the effects on the resulting protein.  

> In simple words:  
> You give the program a gene sequence and the error rates of the viral polymerase, and it tells you:  
> • what % of proteins will stay unchanged  
> • what % will carry amino acid changes  
> • how often premature stop codons appear  
> • what happens in any region of interest you choose  

This can be useful for:
- Understanding how mutation bias affects protein stability  
- Exploring **regions of interest** within a gene  
- Building intuition about error-prone amplification in virology or synthetic biology  

---

## 🧾 Inputs
1. **FASTA file** — mRNA sequence starting at `AUG` and ending at a stop codon (CDS only).  
2. **Bias matrix (CSV)** — a 4×4 table of base substitution probabilities per round.  
   - Header: `A,C,G,U`  
   - Rows: A, C, G, U  
   - Each row must sum to ~1.0.  
3. **Rounds** — how many amplification cycles to simulate (integer).  
4. **Optional region (nt positions)** — a range like `300-600` to check mutation risk in that part of the sequence.  

---

## 📤 Outputs (V1)
- 🧾 **Sequence length** (nt)  
- 🔠 **Protein product (no mutations)** — original translation  
- 🧩 **Codon count** (number of codons in CDS)  
- ✅ **% unchanged protein** (no amino acid changes at all)  
- 📊 **Distribution of k nonsilent mutations** (probability of 0,1,2… codons changing)  
- ⛔ **% premature STOP codon** (protein truncated early)  
- 🎯 **Region of interest** — probability of ≥1 nonsilent mutation in the specified region  

All results are printed in the terminal (future versions will also produce nice HTML/plots 📈).

---

## 🚀 Quickstart

### 1. Clone and set up
```bash
git clone https://github.com/<your-username>/mutantsim.git
cd mutantsim
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e .

---
### 2. Run on the example
python -m mutantsim \
  --fasta examples/example.fasta \
  --bias examples/bias_matrix.csv \
  --rounds 5 \
  --region 4-15
---
### 3. Example Output
Sequence length (nt): 18
Protein (no mutations): MA...
Codons in CDS: 5
% Unchanged protein: 97.500%
% Premature STOP:    0.050%
P(k=0 nonsilent) = 97.500%
P(k=1 nonsilent) = 2.400%
P(k=2 nonsilent) = 0.100%
ROI codons: 1-5
% ≥1 nonsilent in ROI: 2.500%
---

### Respository structure
mutantsim/
  src/mutantsim/      # core code
  examples/           # example FASTA and bias matrix
  tests/              # minimal tests
  pyproject.toml      # package settings
  README.md



