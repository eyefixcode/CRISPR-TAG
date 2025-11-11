import re
from sequence import fetch_region

def _revcomp(s: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(comp)[::-1]

def scan_ngg(chrom: str, center_genomic: int, half: int = 25):
    """
    Find N20-NGG (+strand) and CCN-N20 (-strand) within ±half around center.
    Returns list of dicts: seq20, pam, strand, cut_genomic, distance
    """
    win_start = max(1, center_genomic - half)
    win_end = center_genomic + half
    seq = fetch_region(chrom, win_start, win_end).upper()
    out = []

    # + strand hits: N20 NGG
    for i in range(0, len(seq) - 23 + 1):
        protospacer = seq[i:i+20]
        pam = seq[i+20:i+23]
        if re.match(r"^[ACGT]{20}$", protospacer) and re.match(r"^[ACGT]GG$", pam):
            cut = (win_start + i + 20) - 3  # 3bp upstream of PAM on +
            out.append({
                "seq20": protospacer,
                "pam": pam,
                "strand": "+",
                "cut_genomic": cut,
                "distance": abs(cut - center_genomic),
            })

    # - strand hits: CCN (reverse PAM) then 20nt downstream → reverse-complement
    for i in range(0, len(seq) - 23 + 1):
        pam_plus = seq[i:i+3]
        if re.match(r"^CC[ACGT]$", pam_plus):
            spacer_plus = seq[i+3:i+23]
            if re.match(r"^[ACGT]{20}$", spacer_plus):
                protospacer = _revcomp(spacer_plus)
                pam = _revcomp(pam_plus)
                cut = (win_start + i) + 3  # 3bp upstream of PAM on - (in + coords)
                out.append({
                    "seq20": protospacer,
                    "pam": pam,
                    "strand": "-",
                    "cut_genomic": cut,
                    "distance": abs(cut - center_genomic),
                })

    # keep only guides whose cut site is within the window (paranoid check)
    return [g for g in out if g["distance"] <= half]