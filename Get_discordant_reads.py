import pysam
import argparse
import pandas as pd

# Define and parse command line arguments
parser = argparse.ArgumentParser(description='Extract discordant reads from bam file')
parser.add_argument('-b', '--bamfile', help='Input bam file', required=True)
parser.add_argument('-o', '--outputfile', help='Output bam file', required=True)
parser.add_argument('-h', '--hierfile', help='Family order file', required=True)
args = parser.parse_args()

# Open input bam file
bamfile = pysam.AlignmentFile(args.bamfile, 'rb')

# Read family order file into dataframe
df = pd.read_csv(args.hierfile)
id_list = df['id'].tolist()

# Open output bam file
outputfile = pysam.AlignmentFile(args.outputfile, 'wb', header=bamfile.header)

# Extract discordant reads and write to output bam file
for read in bamfile.fetch():
    if (read.is_proper_pair) == False and ((read.reference_name in id_list) != (read.next_reference_name in id_list)):
        outputfile.write(read)

# Close input and output bam files
bamfile.close()
outputfile.close()
