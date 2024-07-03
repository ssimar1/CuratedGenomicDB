#!/usr/bin/python

# Usage: Usage: python script.py gbk_file multifasta_file output_file
import sys
import warnings
from Bio import SeqIO
from Bio import BiopythonDeprecationWarning
# Suppress BiopythonDeprecationWarning
warnings.simplefilter("ignore", BiopythonDeprecationWarning)

from Bio import pairwise2

from multiprocessing import Pool


def find_matches(ref_id, ref_seq, cds_sequences):
    matches = []
    match_found = False
    for gene_name, cds_seq in cds_sequences.items():
        alignment = pairwise2.align.globalxx(ref_seq, cds_seq, one_alignment_only=True)[0]
        similarity = alignment[2] / max(len(ref_seq), len(cds_seq))
        if similarity >= 0.80:
            if len(ref_seq) != len(cds_seq):
                matches.append((ref_id, gene_name, cds_seq))
                print(f"Found {ref_id}, but diff length")
            else:
                matches.append((ref_id, gene_name, cds_seq))
                print(f"Found {ref_id}")
            match_found = True
    if not match_found:
        print(f"No match found for reference sequence {ref_id}")
    return matches

def main(gbk_file, multifasta_file, output_file):
    print(f"Starting on {gbk_file}")
    # Parse .gbk file
    cds_sequences = {}
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                gene_name = feature.qualifiers.get("gene", [feature.qualifiers.get("locus_tag", ["Unknown"])[0]])[0]
                translation = feature.qualifiers.get("translation", [""])[0]
                cds_sequences[gene_name] = translation

    # Parse multifasta file
    reference_sequences = {}
    for record in SeqIO.parse(multifasta_file, "fasta"):
        reference_sequences[record.id] = str(record.seq)

    # Parallelize the process
    pool = Pool()
    results = [pool.apply_async(find_matches, args=(ref_id, ref_seq, cds_sequences)) for ref_id, ref_seq in reference_sequences.items()]
    pool.close()
    pool.join()

    # Extract results from multiprocessing results
    matches = []
    for result in results:
        matches.extend(result.get())

    # Save matching sequences to a new multifasta file
    with open(output_file, "w") as f:
        print(f"writing {gbk_file} match file")
        for ref_id, gene_name, cds_seq in matches:
            f.write(f">{ref_id}\n{cds_seq}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py gbk_file multifasta_file output_file")
        sys.exit(1)

    gbk_file = sys.argv[1]
    multifasta_file = sys.argv[2]
    output_file = sys.argv[3]



    main(gbk_file, multifasta_file, output_file)