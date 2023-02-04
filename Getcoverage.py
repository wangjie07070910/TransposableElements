import pysam
import numpy as np
import os

directory = "/path/to/directory"
bam_files = [f for f in os.listdir(directory) if f.endswith("rmdup.bam")]

for bam_file in bam_files:
    bam_path = os.path.join(directory, bam_file)
    bam = pysam.AlignmentFile(bam_path, "rb")
    coverage_sum = np.zeros(bam.lengths[0], dtype=float)
    for i in range(bam.nreferences):
        ref = bam.getrname(i)
        coverage = bam.count_coverage(ref)
        coverage_sum += coverage
    mean_coverage = np.mean(coverage_sum)
    print("Mean coverage of {}: {:.2f}".format(bam_file, mean_coverage))
