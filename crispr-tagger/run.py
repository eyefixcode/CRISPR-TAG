# run.py
import sys
from pathlib import Path
from ensembl import get_canonical_transcript, get_codon_sites_genomic
from guides import scan_ngg
from donor import build_donor
from sequence import get_amplicon_window
from primers import design_primers_centered
from io_utils import write_fasta, write_guides_csv, write_primers_csv, write_guides_fasta_for_crispor

# import the callable we just added
from auto_crispor import run_auto_crispor
from subprocess import run, CalledProcessError

def main():
    # --- Inputs ---
    gene = input("Enter Ensembl Gene ID (e.g., ENSG00000008086 for CDKL5): ").strip()
    tag_side = input("Tag side? (5prime/3prime): ").strip().lower()
    if tag_side not in ("5prime","3prime"):
        print("Invalid tag side; defaulting to 3prime")
        tag_side = "3prime"

    out_prefix = f"output_{gene}_{tag_side}"

    # --- Coordinates ---
    tx_id = get_canonical_transcript(gene)  # returns 'ENST...'
    print("Canonical transcript:", tx_id)
    sites = get_codon_sites_genomic(tx_id)  # pass the ID; function returns dict
    print("Genomic sites:\n", sites)

    center = sites["start_codon_genomic"] if tag_side == "5prime" else sites["stop_codon_genomic"]
    print(f"Center site: {center} (strand: {sites['strand']})")

    # --- sgRNAs within ±25 bp ---
    guides = scan_ngg(sites["chrom"], center, half=25)
    print(f"Found {len(guides)} local NGG guides (±25 bp).")

    # --- Donor (200 nt) ---
    donor = build_donor(sites["chrom"], sites["strand"], center, tag_side)
    write_fasta(f"{out_prefix}_donor_200nt.fasta", f"{gene}_{tag_side}_donor", donor)

    # --- Amplicon (±500) + Primers ---
    amplicon_seq, win_start, win_end = get_amplicon_window(sites["chrom"], center, half=500)
    center_idx = center - win_start
    primer_pairs = design_primers_centered(
        amplicon_seq, center_idx,
        product_min=700, product_max=800, num_return=20
    )

    # --- Write local outputs ---
    write_guides_csv(f"{out_prefix}_sgRNAs.csv", guides)
    write_primers_csv(f"{out_prefix}_primers.csv", primer_pairs)
    write_fasta(f"{out_prefix}_amplicon.fasta", f"{sites['chrom']}:{win_start}-{win_end}", amplicon_seq)

    # --- CRISPOR FASTA (preview) ---
    fasta_path = f"{out_prefix}_sgRNAs_for_crispor.fasta"
    write_guides_fasta_for_crispor(fasta_path, sites["chrom"], guides)
    print(f"\nWrote: {fasta_path} (upload/auto-submit to CRISPOR)\n")
    with open(fasta_path) as f:
        print("=== CRISPOR FASTA Preview ===")
        print(f.read())
        print("=============================\n")

    # --- Auto-submit to CRISPOR + download TSV ---
    tsv_path = f"{out_prefix}_crispor_guides.tsv"  # keep unique per gene/side
    try:
        print("Submitting to CRISPOR…")
        run_auto_crispor(fasta_path=fasta_path, out_tsv=tsv_path, headless=True)
        print(f"Downloaded CRISPOR TSV → {tsv_path}")
    except Exception as e:
        print(f"CRISPOR automation failed: {e}")
        sys.exit(2)

    # --- Merge & filter ---
    out_scored = f"{out_prefix}_sgRNAs_scored.csv"
    out_kept   = f"{out_prefix}_sgRNAs_kept.csv"
    try:
        # call merge script via CLI so it stays decoupled
        run([
            sys.executable, "merge_crispor.py",
            f"{out_prefix}_sgRNAs.csv",
            tsv_path,
            out_scored,
            out_kept
        ], check=True)
        print(f"\nWrote: {out_scored}\nWrote: {out_kept}")
    except CalledProcessError as e:
        print("merge_crispor.py failed:", e)
        sys.exit(3)

    print("\nAll done.")

if __name__ == "__main__":
    main()