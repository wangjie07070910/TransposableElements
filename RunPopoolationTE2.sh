#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=8  # 14 physical cores per task
#SBATCH --mem=60G   # 64GB of RAM
#SBATCH --qos=medium
#SBATCH --time=1-16:00:00
#SBATCH --output=%A_%a.Ath_Popte2.stdout
#SBATCH --array=1-8


ml bwa/0.7.17-foss-2018b
ml samtools/1.9-foss-2018b
ref='/scratch-cbe/users/robin.burns/001Popte2/SomaticTEs/TAIR10_Popte2.fasta'
out='/scratch-cbe/users/robin.burns/001Popte2/SomaticTEs/illumina'
raw=/groups/nordborg/raw.data/athaliana/dna/the1001genomes_datasets/001.raw_reads_fastq_gzipped/fastq

#samples=$out/'mypcrfree.txt'
samples=$out/'myacc100bp'
#samples=$out'/myacc76bp'

export l=$SLURM_ARRAY_TASK_ID\p
line=`sed -n $l $samples`
acc=`echo $line | cut -d ' ' -f1`
echo $acc

CONDA_ENVS_PATH='/groups/nordborg/projects/suecica/005scripts/001Software'
CONDA_PKGS_DIRS=/opt/anaconda/pkgs='/groups/nordborg/projects/suecica/005scripts/001Software'
#ml anaconda3/2019.03
#conda create -p /groups/nordborg/projects/suecica/005scripts/001Software/RobinCondaSCRATCH python=2.7 #or 3.6 just run this command once then source it
cd /groups/nordborg/projects/suecica/005scripts/001Software
ml anaconda3/2019.03
source activate RobinCondaSCRATCH3/
cd $out

seqtk trimfq -b 12 -e 12 ${raw}/${acc}.1.fastq.gz > ${out}/${acc}.1.trim.fastq
seqtk trimfq -b 12 -e 12 ${raw}/${acc}.2.fastq.gz > ${out}/${acc}.2.trim.fastq

bwa mem -t 8 -M -U 15 -R '@RG\tID:'$acc'\tSM:'$acc'\tPL:Illumina\tLB:'$acc $ref $out/${acc}.1.trim.fastq $out/${acc}.2.trim.fastq > $out/$acc.sam
#bwa mem -t 8 -M -U 15 -R '@RG\tID:'$acc'\tSM:'$acc'\tPL:Illumina\tLB:'$acc $ref $raw/${acc}.1.fastq.gz $raw/${acc}.2.fastq.gz > $out/$acc.sam

samtools view -@ 8  -bh -t $ref.fai -o $out/$acc.bam $out/$acc.sam
samtools sort -@ 8 -o $out/$acc.sort.bam $out/$acc.bam
samtools index $out/$acc.sort.bam

samtools rmdup $out/$acc.sort.bam $out/$acc.rmdup.bam
samtools index $out/$acc.rmdup.bam
#samtools view -F 256 -h -b $out/$acc.rmdup.bam > $out/$acc.filt.bam  
#samtools index $out/$acc.filt.bam

#count reads mapped afer removing PCR duplicates
#samtools view -c $out/$acc.rmdup.bam > $out/$acc.rmdup_readcount.txt
#Take first 1,000,000 and see what read length is
#samtools view $out/$acc.rmdup.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c > $out/$acc.rmdup_readlength.txt

hier='/groups/nordborg/projects/transposons/Robin/Athaliana/001Popte2/TAIR10_TE_hier.txt'
#ml java/1.8.0_212
popte2='/groups/nordborg/projects/suecica/005scripts/001Software/popte2-v1.10.04.jar'
ml java/1.8.0_212
TMP='/scratch-cbe/users/robin.burns/tmp'
cd $out
java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} ppileup --bam ${acc}.rmdup.bam --map-qual 15 --hier ${hier} --output ${acc}.ppileup.gz
# identify TE signatures
java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} identifySignatures --ppileup ${acc}.ppileup.gz --mode separate --output ${acc}.signatures --min-count 2
# calculate frequencies
java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} frequency --ppileup ${acc}.ppileup.gz --signature ${acc}.signatures --output ${acc}.freqsig
# write TE output 
java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} pairupSignatures --signature ${acc}.freqsig --ref-genome $ref --hier ${hier} --min-distance -200 --max-distance 500 --output ${acc}.teinsertions  

