import pysam
import numpy as np
import os

directory = "/path/to/directory"
bam_files = [f for f in os.listdir(directory) if f.endswith("rmdup.bam")]

for bam_file in bam_files:
    bam_path = os.path.join(directory, bam_file)
    bam = pysam.AlignmentFile(bam_path, "rb")
    max_length = max(bam.lengths[:5])
    coverage_sum = np.zeros((max_length), dtype=float)
    for i in range(5):
        ref = bam.getrname(i)
        coverage = bam.count_coverage(ref)
        coverage_sum[:len(coverage[0])] += np.sum(coverage, axis=0)
    mean_coverage = np.mean(coverage_sum)
    print("Mean coverage of {}: {:.2f}".format(bam_file, mean_coverage))

