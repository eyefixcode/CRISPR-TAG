import requests

ENSEMBL_REST = "https://rest.ensembl.org"

def fetch_region(chrom: str, start: int, end: int) -> str:
    """Fetch 1-based inclusive genomic sequence on the +strand reference."""
    if start < 1:
        start = 1
    url = f"{ENSEMBL_REST}/sequence/region/human/{chrom}:{start}..{end}:1"
    r = requests.get(url, headers={"Accept": "text/plain"}, timeout=30)
    r.raise_for_status()
    return r.text.strip()

def get_amplicon_window(chrom: str, center_genomic: int, half: int = 500):
    """
    Returns (seq, start, end) for a ±half window around center.
    All coords are genomic 1-based inclusive on +strand.
    """
    start = max(1, center_genomic - half)
    end = center_genomic + half
    seq = fetch_region(chrom, start, end)
    return seq.upper(), start, end

