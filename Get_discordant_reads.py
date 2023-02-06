import pysam
import argparse
import csv

def extract_discordant(input_bam, hier_file, output_bam):
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    chrs = ['6244_Chr1', '6244_Chr2', '6244_Chr3', '6244_Chr4', '6244_Chr5']
    id_set = set()
    with open(hier_file, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            id_set.add(row[0])

    outbam = pysam.AlignmentFile(output_bam, "wb", template=bamfile)

    for read in bamfile.fetch():
        if read.is_paired:
            if (read.reference_name in chrs) != (read.next_reference_name in id_set):
                outbam.write(read)

    outbam.close()
    bamfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract discordant reads")
    parser.add_argument("-b", "--bamfile", type=str, required=True, help="Input bam file")
    parser.add_argument("-h", "--hierfile", type=str, required=True, help="Family order csv file")
    parser.add_argument("-o", "--outfile", type=str, required=True, help="Output bam file")
    args = parser.parse_args()

    extract_discordant(args.bamfile, args.hierfile, args.outfile)
