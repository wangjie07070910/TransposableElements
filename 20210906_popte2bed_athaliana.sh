#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=40G   # 72B of RAM
#SBATCH --qos=short
#SBATCH --time=0-06:00:00
#SBATCH --output=%A_%a.popte2bed.stdout
#SBATCH --array=1-15

####
# change the following parameters accordingly 
out='/scratch-cbe/users/robin.burns/001Popte2/SomaticTEs/PCRfree'

index=${SLURM_ARRAY_TASK_ID:-1}
#local_tmpdir='/tmp' # again using local disks
cd $out 
acnum=$(head -n $index $out'/myacc'|tail -n 1) # prefix of file, accession ID from file
accession=${acnum}.teinsertions # accession ID
dist=200 # distance to merge two TEs between individuals
bamfile=${acnum}.sort.bam # bam file used to run popte2 
TAIR10TE=/groups/nordborg/projects/transposons/Robin/Athaliana/001Popte2/TAIR10_TE_hier.txt
outdir=$out #directory to store output files 
date

# cluster TE calls with distance defined above, dist, if the cluster has only one element, there is nothing to do, save in *_cluster_filtered otherwise we will check them *_cluster2check
ml bedtools/2.27.1-foss-2018b
sort -k2,2 -k3,3n -k5,5 ${accession}| awk -v OFS="\t" -v  A=${acnum} '$1=A{print $2,$3,$3+1,$5,$8,$4,$6,$7,$1}' | bedtools cluster -i - -d $dist > ${accession}_cluster
awk '{c[$NF]++}END{for(elem in c){if(c[elem]>1)print elem,c[elem]}}' ${accession}_cluster > ${accession}_cluster2check
awk '{c[$NF]++}END{for(elem in c){if(c[elem]==1)print elem}}' ${accession}_cluster > ${accession}_clustersgood
awk -vCG=${accession}_clustersgood 'BEGIN{while((getline < CG) > 0)c[$1]=1}{if($NF in c)print $0}' ${accession}_cluster > ${accession}_cluster_filtered


ml samtools/1.9-foss-2018b
# look for reads in the vicinity, +- 200 bp. Extract those reads from the bamfile. and the TE reads.
# get all the clusters to check in bed
samtools index $bamfile
awk -vcluster=${accession}_cluster2check 'BEGIN{while((getline < cluster) > 0 )c[$1]=1}{if($NF in c)print $0}' ${accession}_cluster > ${accession}_cluster2check.bed
# get all the reads in all those clusters 
awk -vD=$dist '{if(NR==1){chr=$1;min=$2;max=$3;cl=$NF}else{if(chr==$1 && $NF==cl){if($2<min)min=$2;if($3>max)max=$3}else{print chr":"min-D"-"max+D;chr=$1;min=$2;max=$3;cl=$NF}}}END{print chr":"min-D"-"max+D}' ${accession}_cluster2check.bed |xargs samtools view --threads=10  ${bamfile} > ${acnum}.allreads2check_inregions
awk -vOFS="\t" -vD=$dist '{if(NR==1){chr=$1;min=$2;max=$3;cl=$NF}else{if(chr==$1 && $NF==cl){if($2<min)min=$2;if($3>max)max=$3}else{print chr,min-D,max+D,cl;chr=$1;min=$2;max=$3;cl=$NF}}}END{print chr,min-D,max+D,cl}' ${accession}_cluster2check.bed > ${accession}_cluster2check_cluster600.bed 
# get read names of the reads that overlap a TE
awk '{if($7~/TE/)print $1}' ${acnum}.allreads2check_inregions > ${acnum}.allreads2checkwithTEs_inregions
# get the TE reads with those names
samtools view --threads 10 ${bamfile} | awk -vid=${acnum}.allreads2checkwithTEs_inregions 'BEGIN{while((getline < id) > 0 ){
ID[$1]=1
}}{if($1 in ID){ print $0}}' > ${acnum}.allreads2check_TEs


while read -r line 
do
     cl=$(echo $line |awk '{print $NF}')
     awk -vcluster=${cl} '$NF==cluster'  ${accession}_cluster > ${accession}_cluster.cluster_tmp
     awk -vchr=$(echo $line|awk '{print $1}') -vstart=$(echo $line|awk '{print $2}') -vend=$(echo $line|awk '{print $3}') '$3==chr && $4>=start && $4<=end' ${acnum}.allreads2check_inregions | awk '{if($7~/TE/)print $1}' > ${acnum}.reads2check_inregion
     awk -vid=${acnum}.reads2check_inregion 'BEGIN{while((getline < id) > 0 ){
     ID[$1]=1
     }}{if($1 in ID){ print $0}}' ${acnum}.allreads2check_TEs |sort -k1,1 > ${acnum}.reads2check_TEs
     # check if reads map to only one TE or more, if more, is it from the same family? if not, do we merge calls or not? based on how close and crossmap or not
     # Reads are now sorted by identifier, the reads that are covering the region and their respective alignment on TE were extracted, will keep names of all the secondary hits, and print a table with the start, from the alignment where the first read aligns to the chr, TE names in AT nomenclature. Everything sorted by start, those reads where secondary alignment is another chr are removed, therefore, starts with no reads to TEs are removed. 
     # For the second big block of code I am trying to find the longest path a TE has, paths are 1s, it starts in line 1 with as many 1s as possible and continues, once the path is not there, it doesnt look for that one element anymore until the elements in the initial path array are all cut, then we move to the following cluster. It prints TEid in AT form and number of elements, if the cluster has only one element, it will be removed cause we are asking for at least 2 reads to support the TE. *_clusters.longespath_tmp
        awk '{if(NR==1){read=$1;chrf=0
               if($3~/Chr/){start[read]=$4;
               }else{
                    if($12!~/Chr/){
                         TEatIDcount[$3]++;TEatID[read]=$3
                         mlt[read]=0
                         if($12~/XA:Z/){
                              mlt[read]=1
                              split($12,multi,":");split(multi[3],mtches,";");
                              for(i=1;i<=(length(mtches)-1);i++){
                                   split(mtches[i],tmpte,",");
                                   if(TEatID[read]!=tmpte[1])TEatIDcount[tmpte[1]]++
                                   otherte[read][tmpte[1]]=1
                                   delete tmpte;
                              }
                         }
                    }else{chrf=1}
               }
          }else{
               if($1==read){
                    if($3~/Chr/){start[read]=$4;
                    }else{
                         if($12!~/Chr/){
                              TEatIDcount[$3]++;TEatID[read]=$3
                              mlt[read]=0
                              if($12~/XA:Z/){
                                   mlt[read]=1
                                   split($12,multi,":");split(multi[3],mtches,";");
                                   for(i=1;i<=(length(mtches)-1);i++){
                                        split(mtches[i],tmpte,",");
                                        if(TEatID[read]!=tmpte[1])TEatIDcount[tmpte[1]]++
                                        otherte[read][tmpte[1]]=1
                                        delete tmpte;
                                   }
                              }
                         }else{chrf=1}
                    }
               }else{if(chrf==1 ){delete start[read]}
                    read=$1;chrf=0
                    if($3~/Chr/){start[read]=$4;
                    }else{
                         if($12!~/Chr/){
                              TEatIDcount[$3]++;TEatID[read]=$3
                              mlt[read]=0
                              if($12~/XA:Z/){
                                   mlt[read]=1
                                   split($12,multi,":");split(multi[3],mtches,";");
                                   for(i=1;i<=(length(mtches)-1);i++){
                                        split(mtches[i],tmpte,",");
                                        if(TEatID[read]!=tmpte[1])TEatIDcount[tmpte[1]]++
                                        otherte[read][tmpte[1]]=1
                                        delete tmpte;
                                   }
                              }
                         }else{chrf=1}
                    }
               }
          }
     }END{ printf "%s","start"
                         for(elem in TEatIDcount){
                              if(TEatIDcount[elem]>=2)printf " %s",elem
                         }
                         print ""
                         for(reads in start){
                         printf "%s", start[reads];
                              for(elem in TEatIDcount){
                                   if(TEatIDcount[elem]>=2){
                                        fp=0
                                        if(mlt[reads]==1){
                                             for(tes in otherte[reads]){
                                                  if((tes==elem || elem == TEatID[reads]) && fp==0){printf " %s",1;fp=1}
                                             }
                                        }else{
                                             if(elem == TEatID[reads]){printf " %s",1;fp=1}
                                        }
                                   if(fp==0){printf " %s",0}
                                   }
                              }
                         print "";
                         }
                    }'  ${acnum}.reads2check_TEs |sort -k1,1n|awk '{if(NR==1){print $0}else{for(i=2;i<=NF;i++){s+=$i};if(s>0){print $0};s=0}}'|
     awk '{if(NR==1){for(i=2;i<=NF;i++)TE[i]=$i}else if(NR==2){c++
          for(i=2;i<=NF;i++){if($i==1){path[c][j++];path[c][j]=i;pathc[c][j]++}}
          }else{flagcurrent=0;rmsome=0;
               for(i=2;i<=NF;i++){
                    for(k=1;k<=j;k++){
                         if(path[c][k]==i){
                              if($i==1){ pathc[c][k]++;flagcurrent=1;pathtokeep[c][i]=1;
                              }else{
                              rmsome=1
                              }
                         }
                    }
               }
               if(rmsome==1 && flagcurrent==1){
                    nj=0;
                    for(k=1;k<=j;k++){
                         for(ie in pathtokeep[c]){
                              if(ie==path[c][k]){
                                   tmppath[c][nj++];tmppath[c][nj]=path[c][k];tmppathc[c][nj]=pathc[c][k]
                              }
                         }
                    }
                    delete pathc[c];delete path[c];delete pathtokeep[c]
                    for(k=1;k<=nj;k++){
                         path[c][k];path[c][k]=tmppath[c][k];pathc[c][k]=tmppathc[c][k]
                    }
                    delete tmppath[c]; delete tmppathc[c]
                    j=nj;
               }else if(flagcurrent==0){
                    j=0;c++;
                    for(i=2;i<=NF;i++){if($i==1){path[c][j++];path[c][j]=i;pathc[c][j]++}}
               }
          }
     }END{
          for(i=1;i<=c;i++){
               max[i];maxc=0; max[i]=0;
               for(j=1;j<=length(path[i]);j++){
                    if(maxc<pathc[i][j]){
                         maxc=pathc[i][j];max[i]=TE[path[i][j]];
                    }else if(maxc==pathc[i][j] && max[i]!=TE[path[i][j]]){
                         max[i]=max[i]":"TE[path[i][j]]}
               }
               print i,maxc,max[i]
          }
     }' | awk '$2>1' > ${accession}_clusters.longespath_tmp
     # If there are two equally scored TEs in one cluster, I kept them both separated by colon :. Remove a cluster if the elements in : are not the same family, regardless of same superfamily, if we cant assign it reliably to only one then it might introduce noise. Remove clusters that are == 1 in length, two reads to support insertion, if the same family is found by different TEs, they should be different insertions but as I only see the family, remove to avoid the merging of two different alleles that appear the same it. Trying to avoid allele heterogeneity. As all the lines will be treated similarly, chances are that we will not merge or confound much more than what popte has done already. Append results to filtered file  
     awk -vfnd=${accession}_clusters.longespath_tmp -v TEID=$TAIR10TE 'BEGIN{
          while((getline < TEID) > 0){TE[$1]=$2};
          while((getline < fnd) > 0){
               if($3~/:/){split($3,tes,":");
                    for(elem in tes){
                         TEsame[TE[tes[elem]]]++;
                    };
                    if(length(TEsame)==1){TEfndfam[TE[tes[elem]]][tes[elem]]++};
                    delete TEsame
               }else{TEfndfam[TE[$3]][$3]++;c++
               }
          }
     }{if(($4 in TEfndfam) && length(TEfndfam[$4])==1)print $0}' ${accession}_cluster.cluster_tmp >> ${accession}_cluster_filtered
done < ${accession}_cluster2check_cluster600.bed  

# all the filtered elements are now appended to the *_cluster_filtered file, sort and remove the last column, that contains the cluster number 
sort -k1,1 -k2,2n -k3,3n ${accession}_cluster_filtered | awk 'NF{NF--};1' > ${accession}_cluster_filtered_tmp
#copy filtered file back to lustre
mv ${accession}_cluster_filtered_tmp ${acnum}_cluster_filtered_tmp
#rm ${accession}_clusters.longespath_tmp  ${acnum}_TEstocheck*  ${acnum}*allTEs* ${acnum}*reads2check* 
date
