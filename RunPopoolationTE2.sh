#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=8  # 14 physical cores per task
#SBATCH --mem=84G   # 64GB of RAM
#SBATCH --qos=medium
#SBATCH --time=0-18:01:00
#SBATCH --output=%A_%a.Ath_negcontrol.stdout
# SBATCH --array=1-30


ml bwa/0.7.17-foss-2018b
ml samtools/1.9-foss-2018b
ml spark/2.4.0-hadoop-2.7-java-1.8

refgenome=Popte2_editedgenome.fasta
out=outdir

raw='location/of/fastq/files'

samples=$raw/'myacc'
export l=$SLURM_ARRAY_TASK_ID\p
line=`sed -n $l $samples`
acc=`echo $line | cut -d ' ' -f1`
echo $acc


CONDA_ENVS_PATH='/groups/nordborg/projects/suecica/005scripts/001Software'
CONDA_PKGS_DIRS=/opt/anaconda/pkgs='/groups/nordborg/projects/suecica/005scripts/001Software'
cd /groups/nordborg/projects/suecica/005scripts/001Software
ml anaconda3/2019.03
source activate RobinCondaSCRATCH3/

gunzip $raw/$acc.1.fastq.gz
gunzip $raw/$acc.2.fastq.gz
cd $out
#Remove low complex reads
prinseq-lite.pl -fastq $raw/${acc}.1.fastq -fastq2 $raw/${acc}.2.fastq -out_good ${acc} -out_format 3 -lc_method entropy -lc_threshold 70

bwa mem -t 16 -M -U 15 -R '@RG\tID:'$acc'\tSM:'$acc'\tPL:Illumina\tLB:'$acc $refgenome $out/${acc}_1.fastq $out/${acc}_2.fastq > $out/$acc.sam

samtools view -@ 16  -bh -t $refgenome.fai -o $out/$acc.bam $out/$acc.sam
samtools sort -@ 16 -o $out/$acc.sort.bam $out/$acc.bam
samtools index $out/$acc.sort.bam

#this step isnt needed for pcrfree
ml picard/2.18.27-java-1.8
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        I=$out/$acc.sort.bam \
	      O=$out/$acc.rmdup.bam \
	      M=$out/${acc}_dup_metrics.txt\ 
        REMOVDE_DUPLICATES=true


samtools index $out/$acc.rmdup.bam

#Remove non unique reads

#ILLUM
samtools view -F 256 -q 20 -h -b $out/$acc.rmdup.bam > $out/$acc.filt.bam
samtools index $out/$acc.filt.bam

#PCR
samtools view -F 256 -q 20 -h -b $out/$acc.sort.bam > $out/$acc.filt.bam
samtools index $out/$acc.filt.bam


hier='/location/hier/file'
ml java/1.8.0_212
popte2='/groups/nordborg/projects/suecica/005scripts/001Software/popte2-v1.10.04.jar'
ml java/1.8.0_212
TMP='/scratch-cbe/users/robin.burns/tmp'
mkdir -p $TMP
cd $out
java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} ppileup --bam $out/$acc.filt.bam --map-qual 20 --hier ${hier} --output ${acc}.ppileup.gz
### identify TE signatures
java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} identifySignatures --ppileup ${acc}.ppileup.gz --mode separate --output ${acc}.signatures --min-count 2
## calculate frequencies
java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} frequency --ppileup ${acc}.ppileup.gz --signature ${acc}.signatures --output ${acc}.freqsig
## write TE output 
java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} pairupSignatures --signature ${acc}.freqsig --ref-genome $refgenome --hier ${hier} --min-distance -100 --max-distance 500 --output ${acc}.teinsertions  



#####################################
############Coverage check##############

samtools view -f 3 -h -b $out/${acc}.filt.bam > $out/$acc.filt_properpairs.bam
samtools index $acc.filt_properpairs.bam
samtools view -f 1 -F 256 -F 2 -q 5 -h -b $out/${acc}.filt.bam   > $out/$acc.discF2.bam
samtools index $out2/$acc.discF2.bam

samtools view -h $out/$acc.discF2.bam | awk -F ' ' '$3~"Chr" && $7~"exon" || $3~"exon" && $7~"Chr" || $1 ~ "@" {print $0}' | samtools view  -Sb - > $out/$acc.discF2_chrexon.bam

samtools index $out/$acc.discF2_chrexon.bam


cd $out
ml gcc/8.3.0
python pysam_overlap.py -d $out -m Pop.clustered.exonsmerged.PA.bed -b ${acc}.discF2_chrexon.bam -H exons_400bp.noblast.hier -a ${acc} -p ${acc}.filt_properpairs.bam -c Chr1 
python pysam_overlap.py -d $out -m Pop.clustered.exonsmerged.PA.bed -b ${acc}.discF2_chrexon.bam -H exons_400bp.noblast.hier -a ${acc} -p ${acc}.filt_properpairs.bam -c Chr2
python pysam_overlap.py -d $out -m Pop.clustered.exonsmerged.PA.bed -b ${acc}.discF2_chrexon.bam -H exons_400bp.noblast.hier -a ${acc} -p ${acc}.filt_properpairs.bam -c Chr3
python pysam_overlap.py -d $out -m Pop.clustered.exonsmerged.PA.bed -b ${acc}.discF2_chrexon.bam -H exons_400bp.noblast.hier -a ${acc} -p ${acc}.filt_properpairs.bam -c Chr4
python pysam_overlap.py -d $out -m Pop.clustered.exonsmerged.PA.bed -b ${acc}.discF2_chrexon.bam -H exons_400bp.noblast.hier -a ${acc} -p ${acc}.filt_properpairs.bam -c Chr5


sort -V ${acc}.covfiltChr1.bed > ${acc}.covfiltsortChr1.bed
sort -V ${acc}.covfiltChr2.bed > ${acc}.covfiltsortChr2.bed
sort -V ${acc}.covfiltChr3.bed > ${acc}.covfiltsortChr3.bed
sort -V ${acc}.covfiltChr4.bed > ${acc}.covfiltsortChr4.bed
sort -V ${acc}.covfiltChr5.bed > ${acc}.covfiltsortChr5.bed


cat ${acc}.covfiltsortChr1.bed ${acc}.covfiltsortChr2.bed ${acc}.covfiltsortChr3.bed ${acc}.covfiltsortChr4.bed ${acc}.covfiltsortChr5.bed > ${acc}.covfiltsortChr1to5.bed
