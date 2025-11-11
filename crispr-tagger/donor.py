# donor.py
from sequence import fetch_region

# Fixed 2x Strep tag insert (84 nt)
STREP2_INSERT = (
    "tggagccatccacagttcgaaaaaggtggaggttctggcggtggatcaggtggaagtgcatggtctcaccctcagtttgagaaa"
)
HOM_ARM = 58  # nt
DONOR_LEN = HOM_ARM * 2 + len(STREP2_INSERT)  # 200

def build_donor(chrom: str, strand: int, site: int, tag_side: str) -> str:
    """
    Returns a 200-nt donor oligo in +strand genomic orientation:
      58 bp upstream homology + 84 nt insert + 58 bp downstream homology

    tag_side = "5prime" → insert immediately AFTER the first base of ATG (between base1 and base2)
    tag_side = "3prime" → insert immediately BEFORE the first base of STOP
    'site' is the genomic coordinate (1-based, +strand ref) of:
       - first base of ATG if 5prime
       - first base of STOP if 3prime
    """
    if tag_side not in {"5prime", "3prime"}:
        raise ValueError("tag_side must be '5prime' or '3prime'")

    # For both 5' and 3' we place the insert "before" site on the +strand ref:
    #  - 5': donor sits between site and site+1 (after first base of ATG)
    #  - 3': donor sits before 'site' (first stop base), i.e., between site-1 and site
    # Using homology arms: [site-58 .. site-1] + INSERT + [site .. site+57]
    up = fetch_region(chrom, max(1, site - HOM_ARM), site - 1)
    dn = fetch_region(chrom, site, site + HOM_ARM - 1)

    donor = (up + STREP2_INSERT + dn).upper()
    if len(donor) != DONOR_LEN:
        raise ValueError(f"Donor length is {len(donor)}, expected {DONOR_LEN}")
    # We report donor in reference +strand orientation. For genes on the - strand,
    # this is still fine for synthesis; just document orientation in your workbook.
    return donor