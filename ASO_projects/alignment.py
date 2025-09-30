from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
import subprocess

def clustalo(input_fasta, output_aln):
    cmd = [
        "wsl", "clustalo", 
        "-i", input_fasta, 
        "-o", output_aln, 
        "--outfmt", "clustal", 
        "--force", 
        "-v"
        "--threads=4"
    ]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error:", result.stderr)
        raise RuntimeError("Clustal Omega failed")
    return output_aln

# Example
if __name__ == "__main__":
    aln_file = clustalo("DMD_cdna_isoforms.fasta", "DMD_aligned.aln")
    print("Alignment written to", aln_file)