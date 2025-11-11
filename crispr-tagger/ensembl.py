import requests

ENSEMBL_REST = "https://rest.ensembl.org"

def _headers():
    # Accept is required; a UA helps avoid throttling
    return {
        "Accept": "application/json",
        "User-Agent": "crispr-tagger/0.1 (contact: you@example.com)"
    }

def sanitize_id(stable_id: str) -> str:
    """
    Ensembl REST lookup/id does not accept versioned IDs like ENST...*.9
    This strips any version suffix after a dot.
    """
    return stable_id.split(".", 1)[0]

def get_canonical_transcript(gene_id: str) -> str:
    """Return canonical transcript ID for a given Ensembl Gene ID."""
    gene_id = sanitize_id(gene_id)
    url = f"{ENSEMBL_REST}/lookup/id/{gene_id}?expand=1"
    r = requests.get(url, headers=_headers(), timeout=30)
    r.raise_for_status()
    data = r.json()

    canonical = data.get("canonical_transcript")
    if canonical:
        return canonical

    transcripts = data.get("Transcript", [])
    if not transcripts:
        raise ValueError(f"No transcripts found for gene {gene_id}.")

    transcripts.sort(
        key=lambda t: (t.get("Translation", {}).get("length", 0), t.get("length", 0)),
        reverse=True,
    )
    return transcripts[0]["id"]

def get_sites(transcript_id: str) -> dict:
    """
    Return chrom, strand, and CDS start/end in transcript coordinates.
    (We’ll map to genomic positions in the next step.)
    """
    transcript_id = sanitize_id(transcript_id)
    url = f"{ENSEMBL_REST}/lookup/id/{transcript_id}?expand=1"
    r = requests.get(url, headers=_headers(), timeout=30)
    r.raise_for_status()
    tx = r.json()

    chrom = tx["seq_region_name"]
    strand = tx["strand"]

    if "Translation" not in tx:
        raise ValueError(f"Transcript {transcript_id} has no coding sequence (non-coding).")

    t = tx["Translation"]
    cds_start = t["start"]  # transcript coords (1-based)
    cds_end = t["end"]      # transcript coords (1-based)

    return {
        "chrom": chrom,
        "strand": strand,
        "cds_start_tx": cds_start,
        "cds_end_tx": cds_end,
        "transcript_id": transcript_id,
        "exons": tx.get("Exon", []),  # we’ll use these to map to genomic next
    }
    
# ------- Add below your existing functions in ensembl.py -------

def _order_exons_for_transcript(exons, strand):
    """
    Return exons ordered in transcript (5'→3') direction.
    For +1 strand: smallest genomic start → largest.
    For -1 strand: largest genomic start → smallest.
    """
    exons_sorted = sorted(exons, key=lambda e: e["start"])
    if strand == -1:
        exons_sorted = list(reversed(exons_sorted))
    return exons_sorted

def _build_tx_segments(exons, strand):
    """
    Build a list of segments mapping transcript positions to genomic.
    Each item: (tx_start, tx_end, g_start, g_end, strand)
    """
    ordered = _order_exons_for_transcript(exons, strand)
    segments = []
    tx_pos = 1  # transcript coordinates are 1-based

    for e in ordered:
        length = e["end"] - e["start"] + 1
        # In transcript space, this exon occupies [tx_pos, tx_pos+length-1]
        segments.append((tx_pos, tx_pos + length - 1, e["start"], e["end"]))
        tx_pos += length

    return segments

def _map_tx_to_genomic(tx_pos, segments, strand):
    """
    Map a single transcript position (1-based) to a genomic coordinate (1-based).
    """
    for tstart, tend, gstart, gend in segments:
        if tstart <= tx_pos <= tend:
            offset = tx_pos - tstart
            if strand == 1:
                return gstart + offset
            else:
                # transcript increases while genomic decreases on - strand
                return gend - offset
    raise ValueError(f"Transcript position {tx_pos} is outside exons.")

def get_codon_sites_genomic(transcript_id: str) -> dict:
    """
    Compute genomic coordinates (1-based, on reference +strand coordinates) of:
      - start_codon_genomic: first base of the ATG
      - stop_codon_genomic: first base of the stop codon
    Uses Translation.start/end which are GENOMIC coordinates.
    """
    info = get_sites(transcript_id)  # returns chrom, strand, cds_start_tx, cds_end_tx, exons, transcript_id
    chrom = info["chrom"]
    strand = info["strand"]

    # Re-fetch transcript to read Translation.start/end as GENOMIC coordinates
    # (You could also return them from get_sites if you prefer.)
    import requests
    ENSEMBL_REST = "https://rest.ensembl.org"
    tid = sanitize_id(transcript_id)
    url = f"{ENSEMBL_REST}/lookup/id/{tid}?expand=1"
    r = requests.get(url, headers={"Accept": "application/json"}, timeout=30)
    r.raise_for_status()
    tx = r.json()

    t = tx["Translation"]
    cds_start_gen = t["start"]  # GENOMIC
    cds_end_gen   = t["end"]    # GENOMIC

    if strand == 1:
        start_codon_genomic = cds_start_gen               # first base of ATG
        stop_codon_genomic  = cds_end_gen - 2             # first base of stop (3 bases end the CDS)
    else:
        # On the negative strand, the start codon is at the high end of the CDS.
        start_codon_genomic = cds_end_gen - 2             # first base of ATG in +strand coords
        stop_codon_genomic  = cds_start_gen               # first base of stop in +strand coords

    return {
        "chrom": chrom,
        "strand": strand,
        "transcript_id": tid,
        "start_codon_genomic": start_codon_genomic,
        "stop_codon_genomic": stop_codon_genomic,
    }