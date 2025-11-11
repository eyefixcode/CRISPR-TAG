# merge_crispor.py
import sys
import pandas as pd

SEQ_CANDIDATES = [
    "seq20", "seq", "Guide Sequence", "Guide sequence", "guideSeq",
    "targetSeq", "Target sequence", "Spacer", "Spacer sequence"
]
EFF_CANDIDATES = [
    "efficiency", "Doench 2016", "Doench2016", "Azimuth", "Azimuth 2016",
    "CFD score", "MIT Spec.", "Efficiency Score"
]
OFF01_CANDIDATES = [
    "off_le1mm", "Off-targets (0-1 mismatches)", "Offtargets 0-1MM",
    "0-1 MM Offtargets", "Offtargets_<=1mm", "Off target 0-1mm"
]

def pick_first(cols, candidates):
    for name in candidates:
        if name in cols:
            return name
    return None

def merge_crispor(guides_csv, crispor_tsv, out_scored, out_kept):
    g = pd.read_csv(guides_csv)                      # our local guides
    c = pd.read_csv(crispor_tsv, sep="\t")          # CRISPOR TSV

    
    # --- debug: inspect what we have ---
    print("\n=== Local guides columns ===")
    print(g.columns.tolist())
    print(g.head(3))

    print("\n=== CRISPOR TSV columns ===")
    print(c.columns.tolist())
    print(c.head(3))
    # --- identify columns dynamically ---
    g_seq = pick_first(g.columns, ["seq20", "seq", "sequence", "Spacer"])
    if g_seq is None:
        raise KeyError(f"Could not find a guide-sequence column in {guides_csv}. Got: {list(g.columns)}")

    c_seq = pick_first(c.columns, SEQ_CANDIDATES)
    if c_seq is None:
        # quick debug help
        raise KeyError(f"Could not find a guide-sequence column in {crispor_tsv}. Got: {list(c.columns)}")

    c_eff = pick_first(c.columns, EFF_CANDIDATES)   # may be None
    c_off = pick_first(c.columns, OFF01_CANDIDATES) # may be None

    # --- normalize sequences for a robust join ---
    g["_seq20_norm"] = g[g_seq].astype(str).str.upper().str.replace(r"\s+", "", regex=True)
    c["_seq20_norm"] = c[c_seq].astype(str).str.upper().str.replace(r"\s+", "", regex=True)

    merged = g.merge(c, on="_seq20_norm", how="left", suffixes=("", "_crispor"))

    # --- create unified columns ---
    merged["seq20"] = merged[g_seq]
    if c_eff:
        merged["efficiency"] = pd.to_numeric(merged[c_eff], errors="coerce")
    else:
        merged["efficiency"] = pd.NA
    if c_off:
        # make it integer-like if possible
        merged["off_le1mm"] = pd.to_numeric(merged[c_off], errors="coerce")
    else:
        merged["off_le1mm"] = pd.NA

    # flag the “self-hit counted as off-target” heuristic: off_le1mm == 1
    merged["flag_selfhit"] = merged["off_le1mm"].fillna(-1).eq(1)
    merged.loc[merged["flag_selfhit"], "off_le1mm"] = 0  # treat as zero for filtering

    # --- filtering rules ---
    eff_ok = merged["efficiency"].isna() | (merged["efficiency"] > 10)
    off_ok = merged["off_le1mm"].isna() | (merged["off_le1mm"] == 0) | merged["flag_selfhit"]
    kept = merged[eff_ok & off_ok].copy()

    # tidy output
    cols_order = [col for col in [
        "seq20", "pam", "strand", "cut_genomic", "distance",
        "efficiency", "off_le1mm", "flag_selfhit"
    ] if col in merged.columns]
    merged_out = merged[cols_order].copy()
    kept_out   = kept[cols_order].copy()

    merged_out.to_csv(out_scored, index=False)
    kept_out.to_csv(out_kept, index=False)
    print(f"Scored: {len(merged_out)}, kept: {len(kept_out)}")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python merge_crispor.py <guides_csv> <crispor_tsv> <out_scored> <out_kept>")
        sys.exit(1)
    merge_crispor(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])