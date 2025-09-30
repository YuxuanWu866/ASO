import requests
import pandas as pd
from Bio import Entrez,SeqIO
Entrez.email = "yuw220@ucsd.edu"

def get_transcripts_f_ensembl(gene_name , species ="human"):
    
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?expand=1"
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, headers=headers)
    r.raise_for_status()
    data = r.json()

    transcripts =[]
    for t in data.get("Transcript", []):
        transcripts.append({
            "id":t.get("id"),
            "biotype":t.get("biotype"),
            "length":t.get("length"),
            "tags": t.get("tags",[])
        })
    
    return pd.DataFrame(transcripts)

def co_t(df):
    ##select canonical transcript
    mane = df[df["tags"].apply(lambda x:"MANE_Select" in x if isinstance(x,list) else False)]
    if not mane.empty:
        return mane.iloc[0]["id"]
    
    canonical =df[df["tags"].apply(lambda x:"Ensembl_canonical" in x if isinstance(x,list) else False)]
    if not canonical.empty:
        return canonical.iloc[0]["id"]
    
    protein_coding = df[df["biotype"]=="protein_coding"]
    if not protein_coding.empty:
        return protein_coding.sort_values("length",ascending =False).iloc[0]["id"]
    
    return None

def mrna_ensembl(transcript_id,seq_type="cdna", out_fasta =None):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{transcript_id}?type={seq_type}"
    headers = {"Content-Type": "text/x-fasta"}

    r = requests.get(server + ext, headers=headers)
    if not r.ok:
        return None
    
    fasta_seq = r.text

    if out_fasta:
        with open(out_fasta, "a") as f:
            f.write(fasta_seq)

    return fasta_seq

def download(df,gene_name,seq_type="cdna"):
    out_fasta = f"{gene_name}_{seq_type}_isoforms.fasta"
    open(out_fasta, "w").close()  # clear
    
    for tid in df["id"]:
        seq = mrna_ensembl(tid, seq_type=seq_type, out_fasta=out_fasta)
        if seq:
            print(f"downloaded: {tid}")
        else:
            print(f"not found: {tid}")
    return out_fasta

def generate_unique_kmers(seq, k=20):
   #unique kmer
    kmers = set()
    for i in range(len(seq) - k + 1):
        kmers.add(seq[i:i+k])
    return list(kmers)


if __name__ == "__main__":  
    gene_name = "DMD"

    df = get_transcripts_f_ensembl(gene_name)
    print("all transcripts：")
    print(df.head())

    canonical_id = co_t(df)
    print("\nCanonical transcript:", canonical_id)

    if canonical_id:
        mrna_seq = mrna_ensembl(canonical_id)
        if mrna_seq:
            kmers = generate_unique_kmers(mrna_seq, k=20)
            print(f"generated {len(kmers)} unique 20-mers")
            print("ep:", kmers[:5])
        else:
            print("not found")

    fasta_file = download(df, gene_name, seq_type="cdna")
    print(f"\nAll save to{fasta_file}")