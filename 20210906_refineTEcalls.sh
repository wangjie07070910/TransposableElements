#!/usr/bin/env bash                                                                               
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10  # 14 physical cores per task
#SBATCH --mem=72G   # 72B of RAM
#SBATCH --qos=medium
#SBATCH --time=1-16:00:00
#SBATCH --output=%A_%a.refine.stdout
#SBATCH --array=1-30
# This will be run once for a single process

/bin/hostname
# These commands will be run by each process
# -l flag prepends the task number to each line of output

ml  samtools/0.1.20-foss-2018b
# genomecov file check per accession if there is or there is no match with the TE, at least 2 read pairs
out='/scratch-cbe/users/robin.burns/001Popte2/SomaticTEs/PCRfree/temp'
cd $out
index=${SLURM_ARRAY_TASK_ID:-1}



acc=$(head -n $index $out'/myacc'  | tail -n 1)

TAIR10TE=/groups/nordborg/projects/transposons/Robin/Athaliana/001Popte2/TAIR10_TE_hier.txt
file=2021_AthcheckTES_somatic_merged.genomecov.nohead #delete header of genomecov file manually
ff=${acc}.TEsrefined1
actmp=${acc}.TEann.tmp
bamfile=${acc}.filt.bam


while read -r line;
  do
     TE=$(echo $line| awk '{print $4}') #which TE was found
     awk -v te=$TE '{if($2==te)print $1}' $TAIR10TE > $actmp #TAIR10 of all the TEs with that name
     echo $line | \
       awk '{print $1":"$2-400"-"$3+400}' | # get reads aligned in the region where TE is annotated +- 400 bp
       xargs samtools view $bamfile |    
       awk '{
          el[$7]++; 
          ele[$7,el[$7]]=$4
        } END {
          for (elem in el) {
               for(i = 1; i <= el[elem]; i++) {
                    printf "%s\t",ele[elem,i];
               };
               print el[elem],elem;
          }
        }' | # print per line position where aln reads start, number of lines with the second read mapped to that chr, which chr
       awk -v teannot=$actmp -v Lchr=$(echo "$line"| awk '{print $1}') -v Lstart=$(echo "$line"| awk '{print $2-400}') -v Lend=$(echo "$line"| awk '{print $3+400}')   \
        'BEGIN{
             while ((getline < teannot) > 0){
                  TE[$1]=1};end=0;start=400000000;
          }{
               if($NF in TE){
                    if($1<start){
                         start=$1
                    };
                    if($(NF-2)>end){
                         end=$(NF-2)
                    }
               }
          }END{
               if(start==400000000 && end==0){
                    print Lchr":"Lstart"-"Lend
               }else{
                    print Lchr":"start"-"end+100}
          }'|  # extend the region for retrieving reads using all the different chr in the second read or keep previous. If the number of reads is less than 2, ie 1 or less then NA, 1 supporting.  
        xargs samtools view $bamfile | 
        awk -v teannot=$actmp -v of=$ff -v OFS="\t" -v accession=$acc\
        'BEGIN{
             while ((getline < teannot) > 0){
                  TE[$1]=1
             }
             rightTE=0;
             tot=0;
             tot1=0;
             chrs=0;
        }{
             el[$7]++;
             tot++;
        }END{
             if(tot<1){
                  print "NA" >> of ;
             }else{
                  for (elem in el){
                      if(elem in TE){
                           rightTE+=el[elem];
                      }
                      # if(elem == "="){
                      #     print el[elem]/2,(el[elem]/2)/tot1,elem;
                      # }else{
                      #     print el[elem],el[elem]/tot1,elem;
                      # }
                      if(elem~/Asue_/){
                           chrs+=el[elem]
                      }
                 }
                # print rightTE/tot1;
                # print chrs/tot1;
                if(rightTE>=1 ){
                  if(chrs>rightTE){
                    print "crossmap" >> of
                  }else{
                    print 1 >> of
                  }
                }else{
                     print 0 >> of
                }
             }
          }' # count number of aln with second reads in a right TE chr, same chr, other chr. If number of right TE is bigger than other chr proper signal.
     rm $actmp
done < $file

# copy back result
#cp $ff /lustre/scratch/users/mayela.soto/popoolationTE2/all/
