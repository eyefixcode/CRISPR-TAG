# io_utils.py
import csv

def write_fasta(path: str, header: str, seq: str):
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        # wrap at 70 chars for readability
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i+70] + "\n")

def write_guides_csv(path: str, guides: list[dict]):
    if not guides:
        # create an empty file with headers
        headers = ["seq20","pam","strand","cut_genomic","distance"]
        with open(path, "w", newline="") as fh:
            csv.writer(fh).writerow(headers)
        return
    headers = list(guides[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=headers)
        w.writeheader()
        for row in guides:
            w.writerow(row)
            
            
def write_primers_csv(path: str, primer_pairs: list[dict]):
    import csv
    headers = ["left_seq","right_seq","tm_left","tm_right","gc_left","gc_right",
               "product_size","left_start_in_window","right_end_in_window"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=headers)
        w.writeheader()
        for p in primer_pairs:
            w.writerow(p)
            
def write_guides_fasta_for_crispor(path: str, chrom: str, guides: list[dict]):
    # CRISPOR accepts multi-FASTA; we’ll name headers with locus info for easier merging
    with open(path, "w") as fh:
        for i, g in enumerate(guides, 1):
            header = f"{chrom}|cut={g['cut_genomic']}|strand={g['strand']}|pam={g['pam']}"
            fh.write(f">{header}\n{g['seq20']}\n")