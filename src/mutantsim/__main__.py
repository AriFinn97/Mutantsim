# src/mutantsim/__main__.py
"""
Command-line interface for MutantSim (V1).

Usage example:
python -m mutantsim --fasta examples/example.fasta --bias examples/polio_proxy.csv --rounds 5 --region 4-15

--FASTA reader → load your input sequence.

--Region parser → handle --region 300-600.

--Argparse CLI → so you can run your program from the terminal.

--Pretty printing → nice, human-readable output instead of raw Python dicts"""

import argparse
from pathlib import Path
import numpy as np

from genetics import clean_mrna, translate
from model import load_bias_matrix_csv, matrix_power_per_round
from analytics import summarize

def read_fasta(path: str) -> str:
    """
    Minimal FASTA reader: concatenates non-header lines.
    Accepts sequences with T (DNA); they will be converted to U (RNA) by clean_mrna().
    """
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return clean_mrna("".join(seq_parts))

def parse_region(text: str | None):
    """
    Parse --region like '300-600' into a (start_nt, end_nt) tuple.
    Returns None if not provided.
    """
    if not text:
        return None
    try:
        a, b = text.split("-")
        a, b = int(a), int(b)
        if a < 1 or b < 1:
            raise ValueError
        if b < a:
            a, b = b, a  # auto-swap if reversed
        return (a, b)
    except Exception:
        raise SystemExit("ERROR: --region must look like 300-600 (positive integers).")

def main():
    ap = argparse.ArgumentParser(description="MutantSim — protein outcomes under biased polymerase errors (V1)")
    ap.add_argument("--fasta", required=True, help="Path to mRNA FASTA (CDS: AUG...STOP).")
    ap.add_argument("--bias", required=True, help="4x4 CSV with header A,C,G,U (per-cycle substitution matrix).")
    ap.add_argument("--rounds", type=int, default=1, help="Replication cycles R (default: 1).")
    ap.add_argument("--region", type=str, default=None, help="Nucleotide region of interest, e.g. 300-600 (1-based, inclusive).")
    args = ap.parse_args()

    # --- Load inputs ---
    if not Path(args.fasta).exists():
        raise SystemExit(f"ERROR: FASTA not found: {args.fasta}")
    if not Path(args.bias).exists():
        raise SystemExit(f"ERROR: bias CSV not found: {args.bias}")
    if args.rounds < 0:
        raise SystemExit("ERROR: --rounds must be >= 0")

    seq = read_fasta(args.fasta)          # cleaned (upper, T->U, only A/C/G/U)
    prot = translate(seq)                 # original translation (no mutations)
    M = load_bias_matrix_csv(args.bias)   # per-cycle matrix
    TR = matrix_power_per_round(M, args.rounds)  # effective after R cycles

    roi = parse_region(args.region)

    # --- Run analytics ---
    res = summarize(seq, TR, roi_nt=roi)

    # --- Pretty print ---
    print("=== MutantSim (V1) ===")
    print(f"FASTA:    {args.fasta}")
    print(f"Bias CSV: {args.bias}")
    print(f"Rounds:   {args.rounds}")
    print()
    print(f"Sequence length (nt): {len(seq)}")
    print(f"Protein (no mutations): {prot[:60]}{'...' if len(prot) > 60 else ''}")
    print(f"Codons in CDS: {res['n_codons']}")
    print(f"% Unchanged protein: {100 * res['p_unchanged_protein']:.6f}%")
    print(f"% Premature STOP:    {100 * res['p_premature_stop']:.6f}%")

    # Show first few PMF entries (k=0..min(10, n))
    pmf = res["pmf_k_nonsilent"]
    kmax = min(10, len(pmf) - 1)
    for k in range(kmax + 1):
        print(f"P(k={k} nonsilent) = {100 * float(pmf[k]):.6f}%")

    if res["roi"]:
        s, e = res["roi"]["codon_span"]
        print(f"ROI codons: {s}-{e}")
        print(f"% ≥1 nonsilent in ROI: {100 * res['roi']['p_any_nonsilent']:.6f}%")

if __name__ == "__main__":
    main()
