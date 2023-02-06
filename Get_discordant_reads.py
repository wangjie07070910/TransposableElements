import pysam
import argparse
import csv

parser = argparse.ArgumentParser(description='Filter discordant reads based on hierarchy')
parser.add_argument('-b', '--bamfile', type=str, help='Input bam file', required=True)
parser.add_argument('-o', '--outfile', type=str, help='Output bam file', required=True)
parser.add_argument('-c', '--hierfile', type=str, help='Hierarchy file with id column', required=True)

args = parser.parse_args()

bamfile = args.bamfile
outfile = args.outfile
hierfile = args.hierfile

chromosomes = ['6244_Chr1', '6244_Chr2', '6244_Chr3', '6244_Chr4', '6244_Chr5']

with open(hierfile) as file:
    reader = csv.reader(file)
    hierarchy = [row[0] for row in reader]

bam = pysam.AlignmentFile(bamfile, "rb")
out = pysam.AlignmentFile(outfile, "wb", template=bam)

for read in bam.fetch():
    if read.is_paired:
        if read.reference_name in chromosomes and read.next_reference_name in hierarchy:
            out.write(read)
        elif read.next_reference_name in chromosomes and read.reference_name in hierarchy:
            out.write(read)

bam.close()
out.close()
