# primers.py
from typing import List, Dict
import primer3

def design_primers_centered(template_seq: str, center_index: int,
                            product_min: int = 700, product_max: int = 800,
                            num_return: int = 20) -> List[Dict]:
    """
    template_seq: uppercase DNA string for the whole amplicon window (e.g., 1000 bp)
    center_index: 0-based index within template_seq of the integration site
    Returns a list of primer pair dicts with sequences, Tm, GC, positions, size.
    """
    # Force primer3 to span the integration site by setting SEQUENCE_TARGET
    # around the center (short target length works well).
    target_len = 10
    target_start = max(0, center_index - target_len // 2)

    params = {
        "SEQUENCE_ID": "integration_window",
        "SEQUENCE_TEMPLATE": template_seq,
        "SEQUENCE_TARGET": [target_start, target_len],
        "PRIMER_TASK": "generic",
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, product_max]],
        # Reasonable primer constraints (tweak if too few pairs)
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MIN_TM": 57.0,
        "PRIMER_MAX_TM": 63.0,
        "PRIMER_MIN_GC": 30.0,
        "PRIMER_MAX_GC": 70.0,
        "PRIMER_MAX_POLY_X": 4,
        "PRIMER_EXPLAIN_FLAG": 1
    }

    res = primer3.bindings.designPrimers(params, {})
    pairs = []
    count = res.get("PRIMER_PAIR_NUM_RETURNED", 0)

    for i in range(count):
        Lseq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
        Rseq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        # positions are "start,length" (0-based start in template)
        Lpos, Llen = res[f"PRIMER_LEFT_{i}"]
        Rpos, Rlen = res[f"PRIMER_RIGHT_{i}"]
        pairs.append({
            "left_seq": Lseq,
            "right_seq": Rseq,
            "tm_left": round(res[f"PRIMER_LEFT_{i}_TM"], 2),
            "tm_right": round(res[f"PRIMER_RIGHT_{i}_TM"], 2),
            "gc_left": round(res[f"PRIMER_LEFT_{i}_GC_PERCENT"], 2),
            "gc_right": round(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 2),
            "product_size": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
            "left_start_in_window": Lpos,
            "right_end_in_window": Rpos + Rlen - 1
        })
    return pairs