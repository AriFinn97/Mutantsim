# src/mutantsim/model.py
#m^2;two cycles m^R R cycles
import csv
import numpy as np

BASES = "ACGU"
IDX = {b: i for i, b in enumerate(BASES)}

def load_bias_matrix_csv(path: str) -> np.ndarray:
    """
    Load a 4x4 mutation probability matrix from CSV.
    Header (first row): A,C,G,U
    Rows (in order):   A,C,G,U
    Each row must sum to ~1.0
    """
    rows = []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        header = [h.strip().upper() for h in header]
        if header != list(BASES):
            raise ValueError(f"CSV header must be: {','.join(BASES)}")
        for r in reader:
            rows.append([float(x) for x in r])
    M = np.array(rows, dtype=float)
    if M.shape != (4, 4):
        raise ValueError("Matrix must be 4x4")
    #  quick sanity check
    if not np.allclose(M.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError("Each row of the matrix must sum to 1.0")
    return M

def matrix_power_per_round(M: np.ndarray, rounds: int) -> np.ndarray:
    """Effective mutation matrix after R rounds: M^R."""
    if rounds < 0:
        raise ValueError("rounds must be >= 0")
    if rounds == 0:
        return np.eye(4)
    return np.linalg.matrix_power(M, rounds)
