# src/mutantsim/analytics.py
"""
Analytics for MutantSim (V1).

Inputs:
- seq: RNA sequence (string). We assume CDS is present; we start at FIRST AUG.
- TR:  4x4 effective mutation matrix after R rounds (numpy array), rows/cols A,C,G,U.

Outputs:
- n_codons: number of codons in CDS (from first AUG to before first STOP)
- p_unchanged_protein: probability the protein AA sequence is unchanged
- pmf_k_nonsilent: distribution over K = number of nonsilent codons (Poisson-binomial)
- p_premature_stop: probability any codon becomes a STOP before the original STOP
- roi: optional dict with probability of ≥1 nonsilent codon in a nucleotide region
"""

from __future__ import annotations
from typing import Dict, List, Tuple
import numpy as np

from .genetics import clean_mrna, CODON_TO_AA
from .model import BASES, IDX


def codons_in_cds(seq: str) -> List[str]:
    """
    Return codons from the FIRST AUG until (but not including) the first STOP.
    This keeps the frame consistent with translate().
    """
    s = clean_mrna(seq)

    # Find the first start codon (AUG). If not found, there is no CDS.
    start = s.find("AUG")
    if start == -1:
        return []

    out: List[str] = []
    for i in range(start, len(s) - 2, 3):
        codon = s[i:i + 3]
        if len(codon) < 3:
            break
        aa = CODON_TO_AA.get(codon, "?")
        if aa == "*":
            break  # stop codon: do not include
        out.append(codon)
    return out


def codon_prob_to_codon(orig: str, tgt: str, TR: np.ndarray) -> float:
    """
    Probability that 'orig' codon becomes 'tgt' codon after R rounds,
    assuming independence across the 3 nucleotide positions.
    """
    p = 1.0
    for o, t in zip(orig, tgt):
        p *= TR[IDX[o], IDX[t]]
    return float(p)


def prob_same_amino_acid(orig_codon: str, TR: np.ndarray) -> float:
    """
    Probability that 'orig_codon' ends up encoding the SAME amino acid
    (including the case of no change).
    """
    aa = CODON_TO_AA.get(orig_codon)
    if not aa or aa == "*":
        return 0.0
    total = 0.0
    for tgt, tgt_aa in CODON_TO_AA.items():
        if tgt_aa == aa:
            total += codon_prob_to_codon(orig_codon, tgt, TR)
    return float(total)


def prob_becomes_stop(orig_codon: str, TR: np.ndarray) -> float:
    """Probability that 'orig_codon' becomes ANY stop codon after R rounds."""
    total = 0.0
    for tgt, tgt_aa in CODON_TO_AA.items():
        if tgt_aa == "*":
            total += codon_prob_to_codon(orig_codon, tgt, TR)
    return float(total)


def poisson_binomial_pmf(q: List[float]) -> np.ndarray:
    """
    PMF for K = number of nonsilent codons when each codon i changes with prob q[i].
    Dynamic-programming update; independent but non-identical Bernoulli trials.
    """
    n = len(q)
    pmf = np.zeros(n + 1, dtype=float)
    pmf[0] = 1.0
    for qi in q:
        pmf[1:] = pmf[1:] * (1.0 - qi) + pmf[:-1] * qi
        pmf[0] = pmf[0] * (1.0 - qi)
    s = pmf.sum()
    if s > 0:
        pmf /= s
    return pmf


def summarize(seq: str, TR: np.ndarray, roi_nt: Tuple[int, int] | None = None) -> Dict:
    """
    High-level summary under mutation model TR.

    roi_nt: optional (start_nt, end_nt), 1-based inclusive.
    """
    # 1) Build CDS codon list from FIRST AUG to before STOP
    cods = codons_in_cds(seq)
    n = len(cods)

    # 2) Per-codon probabilities
    p_same = [prob_same_amino_acid(c, TR) for c in cods]  # same AA after mutation
    q = [1.0 - p for p in p_same]                         # nonsilent risk

    # 3) Protein unchanged = all codons remain synonymous
    p_unchanged = float(np.prod([1.0 - qi for qi in q]))

    # 4) Distribution over K = number of nonsilent codons
    pmf = poisson_binomial_pmf(q)

    # 5) Premature STOP anywhere before original STOP
    stop_probs = [prob_becomes_stop(c, TR) for c in cods]
    p_premature_stop = 1.0 - float(np.prod([1.0 - s for s in stop_probs]))

    # 6) Region-of-interest (nt → codon indices)
    roi = None
    if roi_nt:
        a_nt, b_nt = roi_nt
        a_nt = max(1, int(a_nt))
        b_nt = max(a_nt, int(b_nt))
        start_c = (a_nt - 1) // 3
        end_c = min(b_nt // 3, n - 1)
        if end_c >= start_c:
            sub = q[start_c:end_c + 1]
            p_roi_any = 1.0 - float(np.prod([1.0 - qi for qi in sub]))
        else:
            p_roi_any = 0.0
        roi = {"codon_span": (start_c, end_c), "p_any_nonsilent": p_roi_any}

    return {
        "n_codons": n,
        "p_unchanged_protein": p_unchanged,
        "pmf_k_nonsilent": pmf,
        "p_premature_stop": p_premature_stop,
        "roi": roi,
    }
