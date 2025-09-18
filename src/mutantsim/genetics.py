# src/mutantsim/genetics.py
from typing import Dict, List

CODON_TO_AA: Dict[str, str] = {
    "UUU":"F","UUC":"F",
    "UUA":"L","UUG":"L","CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "AUU":"I","AUC":"I","AUA":"I",
    "AUG":"M",
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S","AGU":"S","AGC":"S",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "UAU":"Y","UAC":"Y",
    "CAU":"H","CAC":"H",
    "CAA":"Q","CAG":"Q",
    "AAU":"N","AAC":"N",
    "AAA":"K","AAG":"K",
    "GAU":"D","GAC":"D",
    "GAA":"E","GAG":"E",
    "UGU":"C","UGC":"C",
    "UGG":"W",
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
    "UAA":"*","UAG":"*","UGA":"*",
}
# if we have a DNA sequence it changes the T-->U
def clean_mrna(seq: str) -> str:
    s = seq.strip().upper().replace("T", "U")
    return "".join(ch for ch in s if ch in "ACGU")

def translate(seq: str) -> str:
    """Translate from first AUG to first STOP (CDS-only assumption)."""
    s = clean_mrna(seq)
    start = s.find("AUG")
    if start == -1:
        return "no AUG found"
    aa = []
    # stops the program if the sequence is less than three
    for i in range(start, len(s)-2, 3):
        codon = s[i:i+3]
        if len(codon) < 3:
            break
        # translates the codon using the chart and stops at stop codons
        amino = CODON_TO_AA.get(codon, "?")
        if amino == "*":
            break
        aa.append(amino)
    return "".join(aa)
