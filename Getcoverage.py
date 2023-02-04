import pysam
import os

directory = "/path/to/directory"
bam_files = [f for f in os.listdir(directory) if f.endswith("rmdup.bam")]

for bam_file in bam_files:
    bam_path = os.path.join(directory, bam_file)
    bam = pysam.AlignmentFile(bam_path, "rb")
    mean_coverage = sum(bam.count_coverage(bam.getrname(i)) for i in range(bam.nreferences)) / bam.nreferences
    print("Mean coverage of {}: {:.2f}".format(bam_file, mean_coverage))
