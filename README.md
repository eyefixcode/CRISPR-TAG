# 🧬 CRISPR-TAG
**Automated sgRNA and donor design for tagged gene integration**

CRISPR-TAG is a Python-based pipeline for designing CRISPR/Cas9 guides and donor oligos for site-specific tag integration (e.g., 2× Strep-tag) at the 5′ or 3′ end of a gene of interest.  
The tool integrates **Ensembl**, **CRISPOR**, and **Primer-BLAST** into a single workflow—producing both SnapGene-ready sequences and spreadsheet outputs for downstream cloning or analysis.

---

## 🚀 Features
- **Gene-based automation** – Input any Ensembl Gene ID to fetch canonical transcript coordinates and codon positions.  
- **Context-aware guide scanning** – Finds all SpCas9 NGG PAM sites ± 25 bp around the start or stop codon.  
- **Off-target & efficiency integration** – Merges CRISPOR TSV results and filters for 0–1 mm off-targets and Doench 2016 > 10.  
- **Donor oligo design** – Builds a 200 bp donor containing the 2× Strep-tag sequence flanked by 58 bp homology arms.  
- **Primer design** – Uses Primer3 to generate 700–800 bp PCR amplicons centered on the integration site.  
- **Multi-format export** –  
  - CSV tables for guides & primers  
  - FASTA files for donor, amplicon, and CRISPOR upload  
  - SnapGene-compatible `.dna` annotations *(optional)*  

---

## 🧩 Workflow Overview
1. **Input:** Ensembl Gene ID + desired tag side (`5prime` or `3prime`)  
2. **Fetch gene info:** Retrieve canonical transcript and codon genomic sites  
3. **Scan for guides:** Identify 20 nt NGG sites within ± 25 bp of target  
4. **Design donor:** Insert 2× Strep-tag immediately after start codon or before stop codon  
5. **Design primers:** Create amplification primers ± 500 bp around target  
6. **Score guides (optional):** Upload FASTA → [CRISPOR](https://crispor.gi.ucsc.edu/) → download TSV → merge scores  
7. **Output:** CSV + FASTA files ready for lab use or visualization  

---

## 📂 Example Output

output_tp53_3prime_sgRNAs.csv
output_tp53_3prime_sgRNAs_for_crispor.fasta
output_tp53_3prime_donor_200nt.fasta
output_tp53_3prime_primers.csv
output_tp53_3prime_amplicon.fasta
output_tp53_3prime_sgRNAs_scored.csv
output_tp53_3prime_sgRNAs_kept.csv

---

## 🧠 Example Usage
```bash
# Create virtual environment
python -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Run design pipeline for TP53 3′ tag integration
python run.py
```

Then upload:
output_tp53_3prime_sgRNAs_for_crispor.fasta

to CRISPOR￼ and select:
	•	Genome: Homo sapiens – UCSC Dec 2013 (GRCh38/hg38)
	•	PAM: 20 bp-NGG – Sp Cas9

Download the guides_hg38-unknownLoc.tsv file and merge:
python merge_crispor.py

## 🧪 Dependencies
	•	Python ≥ 3.10
	•	requests￼ – for Ensembl API
	•	primer3-py￼ – primer design
	•	pandas￼ – data merging
	•	biopython￼ – FASTA I/O

(see requirements.txt for exact versions)

## 🧬 Example: TP53 (ENSG00000141510)
python run.py

→ Finds sgRNAs ± 25 bp around TP53 stop codon, builds 2× Strep-tag donor, designs primers, and outputs guide FASTA for CRISPOR.

## 📊 Filtering Rules (merge_crispor.py)

Guides are kept if:
	•	off_le1mm == 0 or flag_selfhit == True
	•	efficiency > 10

This ensures only highly specific, efficient sgRNAs are selected for downstream cloning or validation.

## 🧩 Project Structure
crispr-tagger/
│
├── ensembl.py          # Retrieve canonical transcripts & codon coords
├── guides.py           # Scan for local NGG PAMs
├── donor.py            # Build donor sequence with 2× Strep-tag
├── sequence.py         # Fetch ±500 bp amplicon window
├── primers.py          # Primer3 wrapper for primer design
├── io_utils.py         # Write CSV/FASTA outputs
├── merge_crispor.py    # Merge CRISPOR TSV results & apply filters
└── run.py              # Main pipeline entrypoint

## 🧭 Future Add-Ons
	•	Automatic ± 60 bp flanking sequence export for full Doench 2016 scoring
	•	SnapGene .dna annotation export
	•	Support for Cas12a and SaCas9 PAMs
	•	Batch processing for multiple genes

## 👩‍🔬 Example Results

For CDKL5 (ENSG00000008086), CRISPR-TAG generated:
	•	10 candidate sgRNAs ± 25 bp of stop codon
	•	200 nt donor oligo with 2× Strep-tag
	•	18 primer pairs (700–800 bp amplicons)
	•	All 10 guides passed 0–1 mm off-target filters

⸻

## 👨‍💻 Author

James Settles
Data Scientist | Automation Engineer | Public Health Scientist
⸻


## Flow Chart

<img width="2481" height="3577" alt="Future-friendly variant (tiny local API)-2025-11-11-053440" src="https://github.com/user-attachments/assets/c4036c87-93f5-43e1-9bde-81e665451e2c" />
