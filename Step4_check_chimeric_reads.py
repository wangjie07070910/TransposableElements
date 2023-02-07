import argparse
import pandas as pd
import pysam
import numpy as np

# Parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--matrix", required=True, help="Input matrix file")
parser.add_argument("-b", "--bam", required=True, help="Input bam file")
parser.add_argument("-H", "--hierarchy", required=True, help="Input hierarchy file")
parser.add_argument("-o", "--output", required=True, help="Output file")
args = parser.parse_args()

# Load input files into dataframes
matrix_df = pd.read_csv(args.matrix, sep=',', header=None)
hierarchy_df = pd.read_csv(args.hierarchy, sep='\t')

# Set the column names for matrix_df
matrix_df.columns = ['sample_chr', 'start', 'end', 'gene', 'value']

# Create lists to store the overlapped and non-overlapped rows
myoverlap1 = []
myoverlap2 = []

# Open the bam file using pysam
bam_file = pysam.AlignmentFile(args.bam, "rb")

# Loop through the matrix file
for i, row in matrix_df.iterrows():
    # Check if the chromosome of the read matches the chromosome of the current row
    if row['sample_chr'] in bam_file.references:
        # Find all reads that match the chromosome and are within 1kb of the start of the current row
        overlapping_reads = bam_file.fetch(row['sample_chr'], row['start'] - 1000, row['start'] + 1000)
        overlap = False
        # Loop through the overlapping reads
        for read in overlapping_reads:
            # Get the XA or SA tag
            for tag in read.get_tags():
                if tag[0] == 'XA' or tag[0] == 'SA':
                    gene = tag[1].split(',')[0]
                    # If there is a matching gene in the hierarchy file
                    if hierarchy_df[hierarchy_df['id'].str.contains(gene)].shape[0] > 0:
                        hierarchy_gene = hierarchy_df.loc[hierarchy_df['id'] == gene, 'family'].iloc[0]
                        # If the hierarchy gene matches the gene of the current row
                        if hierarchy_gene == row['gene']:
                            overlap = True
                            break
            if overlap:
                break
        # Store the current row in the appropriate list and modify the value of column 5 if necessary
        if overlap:
            if row['value'] in [0, 'NA']:
                row['value'] = 1
            myoverlap1.append(row)
        else:
            myoverlap2.append(row)

# Close the bam file
bam_file.close()

# Write the overlapped and non-overlapped rows to separate files
pd.DataFrame(myoverlap1).to_csv(args.output + "_overlap1.csv", index=False)
pd.DataFrame(myoverlap2).to_csv(args.output + "_overlap2.csv", index=False)

