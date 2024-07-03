# usage: script.py input_file output_file

import pandas as pd
import sys

# Load the mutation data from the provided TSV file
input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the mutation data
df = pd.read_csv(input_file, sep='\t')

# Create a unique identifier for each mutation
df['Mutation_ID'] = df.apply(lambda row: f"{row['Gene']}_{row['Wildtype']}{row['Position']}{row['Mutation']}", axis=1)

# Get a sorted list of unique mutations
unique_mutations = sorted(df['Mutation_ID'].unique())

# Create an empty matrix with samples as rows and unique mutations as columns
samples = sorted(df['Sample'].unique())
matrix = pd.DataFrame(0, index=samples, columns=unique_mutations)

# Fill the matrix
for _, row in df.iterrows():
    sample = row['Sample']
    mutation_id = row['Mutation_ID']
    matrix.at[sample, mutation_id] = 1

# Save the matrix to a TSV file
matrix.to_csv(output_file, sep='\t')

print(f"Binary presence/absence matrix saved to {output_file}")

