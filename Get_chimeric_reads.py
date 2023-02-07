import argparse
import pysam

parser = argparse.ArgumentParser(description="Split reads from a BAM file")
parser.add_argument("-b", "--infile", help="Path to the input BAM file", required=True)
parser.add_argument("-c", "--idfile", help="Path to the id file", required=True)
parser.add_argument("-o", "--outfile", help="Path to the output BAM file", required=True)

args = parser.parse_args()

bam_file = args.infile
out_bam_file = args.outfile
id_file = args.idfile

bam = pysam.AlignmentFile(bam_file, "rb")
out_bam = pysam.AlignmentFile(out_bam_file, "wb", template=bam)

split_reads = []

# read the ids into a set
ids = set()
with open(id_file, "r") as f:
    for line in f:
        ids.add(line.strip().split("\t")[0])

for read in bam.fetch():
    if read.reference_name == "6244_Chr1" or read.reference_name == "6244_Chr2" or read.reference_name == "6244_Chr3" or read.reference_name == "6244_Chr4" or read.reference_name == "6244_Chr5":
        if read.has_tag("SA"):
            sa_tag = read.get_tag("SA").split(",")
            sa_chrom = sa_tag[0]
            if sa_chrom in ids:
                split_reads.append(read)
                out_bam.write(read)
        elif read.has_tag("XA"):
            xa_tag = read.get_tag("XA").split(",")
            for xa in xa_tag:
                xa_chrom = xa.split(";")[0]
                if xa_chrom in ids:
                    split_reads.append(read)
                    out_bam.write(read)

print("Number of split reads: ", len(split_reads))

out_bam.close()
bam.close()
