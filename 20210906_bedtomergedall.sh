#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1  # 14 physical cores per task
#SBATCH --mem=72G   # 72B of RAM
#SBATCH --qos=medium
#SBATCH --time=0-16:00:00
#SBATCH --output=%A_%a.popte2.stdout

# This will be run once for a single process
# These commands will be run by each process

out='/scratch-cbe/users/robin.burns/001Popte2/SomaticTEs/PCRfree/temp'
#raw='/scratch-cbe/users/robin.burns/001fastq/001DNA'
#samples=$raw'/accessions_list.txt'
#ref='/scratch-cbe/users/robin.burns/Asue_genome.HiCGeneticMap.270919.Popte2.fasta'

####
# change the following parameters accordingly 
popte2resdir=${out} # directory where the output files are stored, work dir
cd $popte2resdir
local_tmpdir='/scratch-cbe/users/robin.burns/tmp' # tmp dir for sort 
dist=400 # distance between two TEs to be merged, between individuals
outgenomecov="2021_AthcheckTES_somatic_merged.genomecov" # file with accession ID, 0s and 1s
outbed="20221_AthcheckTEs_somatic_merged.bed" # bed file for TE insertions 
acclistfile=${out}/myacc # file with the accession IDs used (full path)
#
####
mkdir -p ${local_tmpdir}
ml bedtools/2.27.1-foss-2018b
cat *_cluster_filtered_tmp|
sort -k1,1 -k2,2n -k3,3n -T ${local_tmpdir} | sed 's/ /\t/g' |\
bedtools cluster -i - -d $dist | 
awk -v OFS="\t" '{if(NR==1){ 
  c=1;cluster=$10;name[c]=$4;freq[c]=$5;evid[c]=$8;chr=$1;start[c]=$2;accession[c]=$9;fam[c]=$7;c++; starts=0; freqs=""; evids=""; accessions=""; 
}else{ 
  if(cluster==$10){ 
    start[c]=$2;accession[c]=$9;name[c]=$4;freq[c]=$5;evid[c]=$8;start[c]=$2;accession[c]=$9;fam[c]=$7;c++; 
  }else{ 
    for (elem in name){ 
      names[name[elem]]++; 
    } 
    if(length(names) > 1){ 
      i=1; 
      for(el in names){ 
        if(names[el]==1){ 
          for (ii=1;ii<=length(name);ii++){ 
            if(name[ii] == el){print chr,start[ii], start[ii]+1,name[ii],freq[ii],".",fam[ii],evid[ii],accession[ii]} 
          } 
        }else{ 
          z=0; 
          for (ii=1;ii<=length(name);ii++){ 
            if(name[ii] == el){ 
              pos[i,++z]=ii; 
            } 
          } 
          mz=z; 
          for(z=1; z<=mz;z++){ 
            if(z==1){ 
              min[i] = start[pos[i,z]]; 
              max[i] = start[pos[i,z]]; 
              freqs= freq[pos[i,z]]; 
              evids= evid[pos[i,z]]; 
              accessions= accession[pos[i,z]]; 
            }else{ 
              if(start[pos[i,z]] < min[i]){ min[i] = start[pos[i,z]]}; 
              if(start[pos[i,z]] > max[i]){ max[i] =start[pos[i,z]]}; 
              freqs= freqs";"freq[pos[i,z]]; 
              evids= evids";"evid[pos[i,z]]; 
              accessions= accessions";"accession[pos[i,z]]; 
            } 
          } 
          fam[i]=fam[pos[i,z-1]]; 
          starts=min[i]+int(((max[i]-min[i])/2)+0.5); 
          print chr, starts, starts+1, el,freqs,".",fam[i],evids,accessions; 
          i++; 
        } 
      } 
      delete name; delete freq; delete evid; delete start; delete accession; delete fam;  delete names; delete pos; delete min; delete max; c=1;cluster=$10;name[c]=$4;freq[c]=$5;evid[c]=$8;chr=$1;start[c]=$2;accession[c]=$9;fam[c]=$7;c++;starts=0;freqs=""; evids=""; accessions=""; 
    }else{ 
      for (y=1;y<c;y++){ 
        if(y==1){ 
          min[1]=start[y]; 
          max[1]=start[y]; 
          freqs= freq[y]; 
          evids= evid[y]; 
          accessions= accession[y]; 
        }else{ 
          if(start[y]<min[1]){min[1] = start[y]}; if(start[y]>max[1]){max[1]= start[y]}; 
          freqs= freqs";"freq[y]; 
          evids= evids";"evid[y]; 
          accessions= accessions";"accession[y]; 
        } 
      } ; 
      starts=min[1]+int(((max[1]-min[1])/2)+0.5); 
      {print chr,starts, starts+1,name[1],freqs,".",fam[1],evids,accessions;delete name; delete freq; delete evid; delete start; delete accession; delete fam; delete names; delete min; delete max; c=1;cluster=$10;name[c]=$4;freq[c]=$5;evid[c]=$8;chr=$1;start[c]=$2;accession[c]=$9;fam[c]=$7;c++;starts=0;freqs=""; evids=""; accessions="";} 
    } 
  } 
} 
}END{ 
  for (elem in name){ 
    names[name[elem]]++; 
  } 
  if(length(names) > 1){ 
    i=1; 
    for(el in names){ 
      if(names[el]==1){ 
        for (ii=1;ii<=length(name);ii++){ 
          if(name[ii] == el){print chr,start[ii], start[ii]+1,name[ii],freq[ii],".",fam[ii],evid[ii],accession[ii]} 
        } 
      }else{ 
        z=0; 
        for (ii=1;ii<=length(name);ii++){ 
          if(name[ii] == el){ 
            pos[i,++z]=ii; 
          } 
        } 
        mz=z; 
        for(z=1; z<=mz;z++){ 
          if(z==1){ 
            min[i] = start[pos[i,z]]; 
            max[i] = start[pos[i,z]]; 
            freqs= freq[pos[i,z]]; 
            evids= evid[pos[i,z]]; 
            accessions= accession[pos[i,z]]; 
          }else{ 
            if(start[pos[i,z]] < min[i]){ min[i] = start[pos[i,z]]}; 
            if(start[pos[i,z]] > max[i]){ max[i] =start[pos[i,z]]}; 
            freqs= freqs";"freq[pos[i,z]]; 
            evids= evids";"evid[pos[i,z]]; 
            accessions= accessions";"accession[pos[i,z]]; 
          } 
        } 
        fam[i]=fam[pos[i,z-1]]; 
        starts=min[i]+int(((max[i]-min[i])/2)+0.5); 
        print chr, starts, starts+1, el,freqs,".",fam[i],evids,accessions; 
        i++; 
      } 
    }  
    delete name; delete freq; delete evid; delete start; delete accession; delete fam;  delete names; delete pos; delete min; delete max; c=1;cluster=$10;name[c]=$4;freq[c]=$5;evid[c]=$8;chr=$1;start[c]=$2;accession[c]=$9;fam[c]=$7;c++;starts=0;freqs=""; evids=""; accessions=""; 
  }else{ 
    for (y=1;y<c;y++){ 
      if(y==1){ 
        min[1]=start[y]; 
        max[1]=start[y]; 
        freqs= freq[y]; 
        evids= evid[y]; 
        accessions= accession[y]; 
      }else{ 
        if(start[y]<min[1]){min[1] = start[y]}; if(start[y]>max[1]){max[1]= start[y]}; 
        freqs= freqs";"freq[y]; 
        evids= evids";"evid[y]; 
        accessions= accessions";"accession[y]; 
      } 
    } ; 
    starts=min[1]+int(((max[1]-min[1])/2)+0.5); 
    {print chr,starts, starts+1,name[1],freqs,".",fam[1],evids,accessions;delete name; delete freq; delete evid; delete start; delete accession; delete fam; delete names; delete min; delete max; c=1;cluster=$10;name[c]=$4;freq[c]=$5;evid[c]=$8;chr=$1;start[c]=$2;accession[c]=$9;fam[c]=$7;c++;starts=0;freqs=""; evids=""; accessions="";} 
  } 
}' > $outbed
 
 sort -k1,1 -k2,2n -k3,3n -T ${local_tmpdir} $outbed | \
 awk -v OFS="\t" -vACC=${acclistfile} 'BEGIN{
   while((getline<ACC)>0){
     a[$1]=1
     };
     printf "%s %s %s %s ","chr","start","end","TE";
     for (i in a) {
       printf "%s ",i
       };
       print ""
       }{
         printf "%s %s %s %s ",$1,$2,$3,$4;
         split($9,acc,";");
         for(i in a){
           fl=0; 
           for(k in acc){
             if(i==acc[k] && fl==0){
               printf "%s ",1;
               fl=1;
               break
               }
               };
               if(fl==0){
                 printf "%s ",0
                 }
                 };
                 print ""
        }'  >  ${outgenomecov}
        
# genomecov to bed, cluster to check easier for reads, basically will re-do the first step but including more TEs 
awk -vOFS="\t" '{if(NR>1){print $1,$2,$3,$4}}' $outgenomecov > $outbed
