"""input fasta align_file target_id
  output : 
20-mer in target columns 
"""
import sys
import csv
from Bio import SeqIO,AlignIO
from collections import Counter

a_fasta = "DMD_cdna_isoforms.fasta"
aln = "DMD_aligned.aln"
align_fmt = "clustal"
target_id = "ENST00000684292.1"
k=20
out_csv=f"{target_id}_unique_{k}mers.csv"
"""
fasta load
"""
seqs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(a_fasta, "fasta")}
if target_id not in seqs:
    tgt_core = target_id.split(".")[0]
    matches = [k for k in seqs if k.split(".")[0] == tgt_core]
    if matches: target_id = matches[0]
    else:
        print("Error: target id not found")
        print(list(seqs.keys())[:10])
        sys.exit(1)

target_seq = seqs[target_id]
other_seqs = {k:v for k,v in seqs.items() if k != target_id}
print(f"Loaded {len(seqs)} sequences, target: {target_id}, length={len(target_seq)}")

"""
k-mers set
"""
other_kmers = set()
for sid, s in other_seqs.items():
    for i in range(len(s) - k + 1):
        other_kmers.add(s[i:i+k])
print(f"Built other isoforms k-mer set (size={len(other_kmers)})")

# load alignment
alignment = AlignIO.read(aln, align_fmt)
aln_ids = [rec.id for rec in alignment]
if target_id not in aln_ids:
    tgt_core = target_id.split(".")[0]
    match = None
    for aid in aln_ids:
        if aid.split(".")[0] == tgt_core:
            match = aid
            break
    if match is None:
        print("ERROR: target id not found in alignment headers. Some alignment headers:")
        print(aln_ids[:20])
        sys.exit(1)
    else:
        target_aln_id = match
else:
    target_aln_id = target_id

print(f"Using alignment ID: {target_aln_id}")

aln_len = alignment.get_alignment_length()
num_seqs = len(alignment)
gap_sensitive_cols = []
row_index = aln_ids.index(target_aln_id)

for col in range(aln_len):
    col_chars = [alignment[r][col].upper() for r in range(num_seqs)]
    non_gaps = [c for c in col_chars if c != '-']
    if non_gaps:
        common, cnt = Counter(non_gaps).most_common(1)[0]
        _ = cnt / len(non_gaps)
    tgt_char = alignment[row_index][col].upper()
    if tgt_char != '-' and any(alignment[r][col] == '-' for r in range(num_seqs) if r != row_index):
        gap_sensitive_cols.append(col)


aln_to_target_pos = {}
tpos = 0
for col in range(aln_len):
    if alignment[row_index][col] != '-':
        tpos += 1
        aln_to_target_pos[col] = tpos
    else:
        aln_to_target_pos[col] = None

targetpos_to_alncol = {pos: col for col, pos in aln_to_target_pos.items() if pos is not None}

# for all k-mers in target
results = []
for i in range(len(target_seq) - k + 1):
    start = i + 1
    kmer = target_seq[i:i+k]
    unique = (kmer not in other_kmers)
    has_N = 'N' in kmer
    gc = (kmer.count('G') + kmer.count('C')) / k

    aln_cols = [targetpos_to_alncol[pos] for pos in range(start, start+k)]
    junction_spanning = any(col in gap_sensitive_cols for col in aln_cols if col is not None)

    results.append({
        "transcript_id": target_id,
        "start_1based": start,
        "kmer": kmer,
        "unique_in_others": unique,
        "junction_spanning": junction_spanning,
        "gc_frac": round(gc,3),
        "has_N": has_N
    })

# save CSV
with open(out_csv, "w", newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=[
        "transcript_id","start_1based","kmer","unique_in_others","junction_spanning","gc_frac","has_N"
    ])
    writer.writeheader()
    for row in results:
        writer.writerow(row)

print(f"Wrote {len(results)} k-mers to {out_csv}")
total_unique = sum(1 for r in results if r["unique_in_others"])
junction_and_unique = sum(1 for r in results if r["unique_in_others"] and r["junction_spanning"])
print(f"Unique {k}-mers: {total_unique}; unique & junction-spanning: {junction_and_unique}")

print("Examples (unique & junction):")
count = 0
for r in results:
    if r["unique_in_others"] and r["junction_spanning"]:
        print(r["start_1based"], r["kmer"], r["gc_frac"])
        count += 1
        if count >= 20: break
