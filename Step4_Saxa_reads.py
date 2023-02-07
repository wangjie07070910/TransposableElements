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

# Open the bam file using pysam
bam_file = pysam.AlignmentFile(args.bam, "rb")

# Loop through the bam file
for read in bam_file.fetch():
    #print("Processing read")
    #print(read)
    #break
    # Check if the chromosome matches
    if read.reference_name == matrix_df.loc[0, 'sample_chr']:
        #print(read.reference_start)
        #print(matrix_df.loc[2,'start'])
        #print(abs(read.reference_start - matrix_df.loc[4,'start']))
        #break
        # Check if the start of the read is within 1kb of col2 of the matrix file
        #if np.any(np.abs(read.reference_start - matrix_df.loc[matrix_df['sample_chr'] == bam_file.get_reference_name(read.reference_id), 'start'].iloc[0]) <= 1000):
        match = np.abs(read.reference_start - matrix_df.loc[,'start']) <= 2000
        if np.any(match):
            # Get the XA or SA tag
            for tag in read.get_tags():
                print("Found a tag")
                #print(tag)
                #break
                if tag[0] == 'XA' or tag[0] == 'SA':
                    gene = tag[1].split(',')[0]
                    print(gene)
                    # If there is a string that matches column1 of the hierarchy file
                    #if gene in hierarchy_df['id'].values:
                    if hierarchy_df[hierarchy_df['id'].str.contains(gene)].shape[0] > 0:
                        hierarchy_gene = hierarchy_df.loc[hierarchy_df['id'] == gene, 'family'].iloc[0]
                        # If Column 2 of the hierarchy file matches col4 in the matrix file
                        if hierarchy_gene in matrix_df['gene'].values:
                            matrix_index = matrix_df.index[matrix_df['gene'] == hierarchy_gene].tolist()[0]
                            # If column5 of the matrix file is 0 or NA, change this to 1
                            if matrix_df.loc[matrix_index, 'value'] in [0, 'NA']:
                                matrix_df.at[matrix_index, 'value'] = 1

# Close the bam file
bam_file.close()

# Write the output to a file
matrix_df.to_csv(args.output, sep=',', index=False)

