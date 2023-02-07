#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1 # 14 physical cores per task
#SBATCH --mem=44G   # 64GB of RAM
#SBATCH --qos=medium
#SBATCH --time=0-16:00:00
#SBATCH --output=%A_%a.Popte2.stdout
#SBATCH --array=1-16


ml bwa/0.7.17-foss-2018b
ml samtools/1.9-foss-2018b
ml spark/2.4.0-hadoop-2.7-java-1.8

#refgenome='/groups/nordborg/projects/transposons/Robin/001TransposableElements/010For1001g/'
#refgenome='/groups/nordborg/projects/transposons/Robin/001TransposableElements/001Popte2/008FilteredTEs/6244.Popte2_010223.fasta'
refgenome='/scratch-cbe/users/robin.burns/27genomes/Popte2/010For1001g/6244_blastmasked_Popte2.fasta'
#fullgenome='/scratch-cbe/users/robin.burns/27genomes/Popte2/010For1001g/6244.fulllibrary.fasta'
#index=/groups/nordborg/projects/transposons/Robin/001TransposableElements/010For1001g/6244.forPopte2
fq='/scratch-cbe/users/robin.burns/27genomes/Popte2'
out='/scratch-cbe/users/robin.burns/27genomes/Popte2'
raw=/groups/nordborg/raw.data/athaliana/dna/the1001genomes_datasets/001.raw_reads_fastq_gzipped/fastq


samples='/groups/nordborg/projects/transposons/Robin/001TransposableElements/002scripts/002Popte2/my16acc.txt'
acc=$(awk "NR==$SLURM_ARRAY_TASK_ID" $samples)
echo $acc

#acc='6124'
CONDA_ENVS_PATH='/groups/nordborg/projects/suecica/005scripts/001Software'
CONDA_PKGS_DIRS=/opt/anaconda/pkgs='/groups/nordborg/projects/suecica/005scripts/001Software'
cd /groups/nordborg/projects/suecica/005scripts/001Software
ml anaconda3/2019.03
#source activate RobinConda_Oct22/
source activate RobinCondaSCRATCH3/

#gunzip $out/$acc.1.fastq.gz
#gunzip $out/$acc.2.fastq.gz

#ml fastqc/0.11.8-java-1.8
#fastqc $out/$acc.1.fastq


#Remove low complex reads
cd $out
#prinseq-lite.pl -fastq $out/${acc}.1.fastq -fastq2 $out/${acc}.2.fastq -out_bad null -out_good ${acc}_prinseq -out_format 3 -lc_method entropy -lc_threshold 70

#acc=8236
#bwa mem -t 16 -M -R '@RG\tID:'$acc'\tSM:'$acc'\tPL:Illumina\tLB:'$acc $refgenome $fq/${acc}_prinseq_1.fastq $fq/${acc}_prinseq_2.fastq > $out/$acc.sam

#bwa mem -t 16 -M -R '@RG\tID:'$acc'\tSM:'$acc'\tPL:Illumina\tLB:'$acc $refgenome $raw/${acc}.1.fastq.gz $raw/${acc}.2.fastq.gz > $out/$acc.sam


#is bowtie2 better?
#ml bowtie2/2.3.4.2-foss-2018b
#acc=8236
#bowtie2 --no-unal -N 1 --very-sensitive-local -D 20 -R 3 -N 1 -L 20 -p 6 -x $index -1 $fq/${acc}_prinseq_1.fastq -2 $fq/${acc}_prinseq_2.fastq -S $out/$acc.sam

#samtools view -@ 16  -bh -t $refgenome.fai -o $out/$acc.bam $out/$acc.sam
#samtools sort -@ 16 -o $out/$acc.sort.bam $out/$acc.bam
#samtools index $out/$acc.sort.bam

#this step isnt needed for pcrfree
#ml picard/2.18.27-java-1.8
#java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
#        I=$out/$acc.sort.bam \
#	O=$out/$acc.rmdup.bam \
#	M=$out/${acc}_dup_metrics.txt \
#     	REMOVE_DUPLICATES=true
#

#samtools index $out/$acc.rmdup.bam

#Remove non unique reads

cd $out
#mysb='/groups/nordborg/projects/suecica/005scripts/001Software/samblaster'

#samtools view -@ 16 -h ${acc}.bam | ${mysb}/samblaster -e -u "${acc}.umap.fastq" | samtools view -@ 16 -b - > ${acc}_unsort.bam

#yhindex='/groups/nordborg/projects/transposons/Robin/001TransposableElements/001Popte2/008FilteredTEs/6244.Popte2_010223.X15_01_65525S'
#proc=16
#yaha -t $proc -x $yhindex -q "${acc}.umap.fastq" -L 11 -H 2000 -M 15 -osh stdout | ${mysb}/samblaster -s "${acc}.split.sam" > /dev/null

#samtools view -@ $proc -bh -t $refgenome.fai -o ${acc}.split.bam ${acc}.split.sam
#samtools sort -@ $proc -o ${acc}.split.sort.bam ${acc}.split.bam
#samtools index ${acc}.split.sort.bam

#cd $out

#ILLUM
#samtools view -F 256 -q 5 -h -b $out/$acc.rmdup.bam > $out/$acc.filt.bam
#samtools index $out/$acc.filt.bam


#hier='/groups/nordborg/projects/transposons/Robin/001TransposableElements/001Popte2/008FilteredTEs/TAIR10_TE_hierarchy_extended.filt130422.forPopte2.txt'
hier='/groups/nordborg/projects/transposons/Robin/001TransposableElements/001Popte2/008FilteredTEs/Ath_RepBase_TEs_filt010223_extended.txt'
#ml java/1.8.0_212
popte2='/groups/nordborg/projects/suecica/005scripts/001Software/popte2-v1.10.04.jar'
ml java/1.8.0_212
TMP='/scratch-cbe/users/robin.burns/tmp'
mkdir -p $TMP
cd $out
#java -Xmx80g -Djava.io.tmpdir=$TMP -jar ${popte2} ppileup --bam $out/$acc.rmdup.bam --map-qual 5 --hier ${hier} --output ${acc}.ppileup.gz
### identify TE signatures
#java -Xmx80g -Djava.io.tmpdir=$TMP -jar ${popte2} identifySignatures --ppileup ${acc}.ppileup.gz --mode separate --output ${acc}.signatures --min-count 2
### calculate frequencies
#java -Xmx80g -Djava.io.tmpdir=$TMP -jar ${popte2} frequency --ppileup ${acc}.ppileup.gz --signature ${acc}.signatures --output ${acc}.freqsig
### write TE output 
#java -Xmx80g -Djava.io.tmpdir=$TMP -jar ${popte2} pairupSignatures --signature ${acc}.freqsig --ref-genome $refgenome --hier ${hier} --min-distance -100 --max-distance 500 --output ${acc}.teinsertions  



#####################################
############Get discordant reads##############

#samtools view -f 3 -h -b $out/${acc}.rmdup.bam > $out/$acc.filt_properpairs.bam
#samtools index $acc.filt_properpairs.bam

scripts='/groups/nordborg/projects/transposons/Robin/001TransposableElements/002scripts/002Popte2'
cd $out
#python ${scripts}/Get_discordant_reads.py -b ${acc}.rmdup.bam -o ${acc}.discordant_TEreads.bam -c Ath_RepBase_TEs_filt010223_extended.txt
#samtools index ${acc}.discordant_TEreads.bam

#Change saxa to chimeric todo
#python ${scripts}/get_saxa_reads.py -b $acc.rmdup.bam -o ${acc}.saxa_TEreads.bam -c Ath_RepBase_TEs_filt010223_extended.txt
#samtools index ${acc}.saxa_TEreads.bam
###################


#############Run Popte2 Step1 and Step2#############


#For Step3
cd $out
scripts='/groups/nordborg/projects/transposons/Robin/001TransposableElements/002scripts/002Popte2'
ml gcc/8.3.0
hier='/groups/nordborg/projects/transposons/Robin/001TransposableElements/001Popte2/008FilteredTEs/Ath_RepBase_TEs_filt010223_extended.txt'
#python ${scripts}/Step3_covcheck_pysam.py -d $out -m Pop.clustered.merged.PA.bed -b ${acc}.discordant_TEreads.bam -H ${hier} -a ${acc} -p ${acc}.filt_properpairs.bam -c 6244_Chr1 
#python ${scripts}/Step3_covcheck_pysam.py -d $out -m Pop.clustered.merged.PA.bed -b ${acc}.discordant_TEreads.bam -H ${hier} -a ${acc} -p ${acc}.filt_properpairs.bam -c 6244_Chr2
#python ${scripts}/Step3_covcheck_pysam.py -d $out -m Pop.clustered.merged.PA.bed -b ${acc}.discordant_TEreads.bam -H ${hier} -a ${acc} -p ${acc}.filt_properpairs.bam -c 6244_Chr3
#python ${scripts}/Step3_covcheck_pysam.py -d $out -m Pop.clustered.merged.PA.bed -b ${acc}.discordant_TEreads.bam -H ${hier} -a ${acc} -p ${acc}.filt_properpairs.bam -c 6244_Chr4
#python ${scripts}/Step3_covcheck_pysam.py -d $out -m Pop.clustered.merged.PA.bed -b ${acc}.discordant_TEreads.bam -H ${hier} -a ${acc} -p ${acc}.filt_properpairs.bam -c 6244_Chr5

cd $out
#sort -V ${acc}.covfilt*Chr1.bed > ${acc}.covfiltsortChr1.bed
#sort -V ${acc}.covfilt*Chr2.bed > ${acc}.covfiltsortChr2.bed
#sort -V ${acc}.covfilt*Chr3.bed > ${acc}.covfiltsortChr3.bed
#sort -V ${acc}.covfilt*Chr4.bed > ${acc}.covfiltsortChr4.bed
#sort -V ${acc}.covfilt*Chr5.bed > ${acc}.covfiltsortChr5.bed


#cat ${acc}.covfiltsortChr1.bed ${acc}.covfiltsortChr2.bed ${acc}.covfiltsortChr3.bed ${acc}.covfiltsortChr4.bed ${acc}.covfiltsortChr5.bed > ${acc}.covfiltsortChr1to5.bed

python ${scripts}/Step4_Saxa_reads.py -H ${hier} -m ${acc}.covfiltsortChr1to5.bed -b ${acc}.saxa_TEreads.bam -o ${acc}.covfiltsortChr1to5_doublecheck.bed 
