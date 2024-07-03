import os
import sys
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MuscleCommandline
import pandas as pd
import time

def run_muscle(input_fasta, output_fasta):
    muscle_executable = "/usr/local/bin/muscle"  # Replace with the full path to MUSCLE if different
    muscle_cline = MuscleCommandline(cmd=muscle_executable, input=input_fasta, out=output_fasta)
    stdout, stderr = muscle_cline()
    print("MUSCLE output:", stdout)
    print("MUSCLE error:", stderr)

def create_combined_fasta(reference_seq, target_seqs, combined_fasta):
    sequences = [reference_seq]
    sequences.extend(target_seqs)
    SeqIO.write(sequences, combined_fasta, "fasta")

def find_differences(msa_file, gene_id, output_tsv):
    alignment = AlignIO.read(msa_file, "fasta")
    differences = []
    
    ref_seq = None
    sample_seqs = []

    for record in alignment:
        if record.id == gene_id:
            ref_seq = record.seq
        else:
            sample_seqs.append(record)

    if ref_seq is None:
        print(f"Reference sequence for {gene_id} not found in alignment.")
        return []

    for record in sample_seqs:
        sample_id = record.id.split('|')[-1]
        match_seq = record.seq

        for pos, (ref_base, match_base) in enumerate(zip(ref_seq, match_seq)):
            if ref_base != match_base:
                if ref_base == '-':
                    differences.append((gene_id, pos + 1, '-', match_base, "Insertion", sample_id))
                elif match_base == '-':
                    differences.append((gene_id, pos + 1, ref_base, '-', "Deletion", sample_id))
                else:
                    differences.append((gene_id, pos + 1, ref_base, match_base, "Substitution", sample_id))

    return differences

def main(reference_fasta, matches_fasta_list, output_tsv):
    start_time = time.time()

    # Read reference sequences
    reference_sequences = {record.id: record for record in SeqIO.parse(reference_fasta, "fasta")}

    all_differences = []

    # Read match sequences from each file
    with open(matches_fasta_list, "r") as f:
        sample_files = [line.strip() for line in f]

    for gene_id, ref_seq in reference_sequences.items():
        print(f"Processing gene: {gene_id}")
        target_seqs = []

        for sample_file in sample_files:
            sample_id = os.path.basename(sample_file).split('.')[0]
            found = False
            for record in SeqIO.parse(sample_file, "fasta"):
                if record.id == gene_id:
                    record.id = f"{record.id}|{sample_id}"
                    record.description = ""
                    target_seqs.append(record)
                    found = True

            if not found:
                print(f"No target sequences found for gene {gene_id} in sample {sample_id}. Adding MISSING entries.")
                differences = [(gene_id, 'MISSING', 'MISSING', 'MISSING', 'MISSING', sample_id)]
                all_differences.extend(differences)

        if not target_seqs:
            print(f"No target sequences found for gene {gene_id} in any samples. Skipping MSA.")
            continue

        combined_fasta = f"combined_{gene_id}.fasta"
        msa_output_fasta = f"msa_output_{gene_id}.fasta"

        create_combined_fasta(ref_seq, target_seqs, combined_fasta)
        run_muscle(combined_fasta, msa_output_fasta)

        differences = find_differences(msa_output_fasta, gene_id, output_tsv)
        all_differences.extend(differences)

    df = pd.DataFrame(all_differences, columns=["Gene", "Position", "Wildtype", "Mutation", "Type", "Sample"])
    df.to_csv(output_tsv, sep='\t', index=False)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total elapsed time: {elapsed_time / 60:.2f} minutes")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py reference_fasta matches_fasta_list output_tsv")
        sys.exit(1)

    reference_fasta = sys.argv[1]
    matches_fasta_list = sys.argv[2]
    output_tsv = sys.argv[3]

    main(reference_fasta, matches_fasta_list, output_tsv)

