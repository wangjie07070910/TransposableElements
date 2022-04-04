#!/usr/bin/env python
import numpy as np
from optparse import OptionParser                                               

# pipeline params
parser = OptionParser()
parser.add_option("-p", "--pairs", dest='pairfile', help='pair read file', default="")
parser.add_option("-b", "--bam", dest="bamfile", help="split bam file", default="")
parser.add_option("-d", "--dir", dest="dir", help="directory of combined vcf file", default="")
parser.add_option("-m", "--matrix", dest="matrix", help="te matrix file", default="")
parser.add_option("-H", "--hier", dest="hier", help="hierarchy file", default="")
parser.add_option("-a", "--acc", dest="acc", help="accession", default="")
parser.add_option("-c", "--chr", dest="chr", help="chromosome", default="")
(options, args) = parser.parse_args()                                            


import os
import pandas as pd
import pysam
import bioframe as bf
from collections import defaultdict
import pybedtools 
splitbam=pysam.AlignmentFile("%s/%s"%(options.dir, options.bamfile), "rb")
pairbam=pysam.AlignmentFile("%s/%s"%(options.dir, options.pairfile), 'rb')

myhier = defaultdict(list)
for line in open("%s/%s"%(options.dir, options.hier), 'r'):
	if 'id' in line:
		pass
	else:
		line=line.strip().split('\t')	
		key,value=line[2], line[0]
		if key not in myhier:
			myhier[key]=[value]
		else:
			myhier[key].append(value)


myreads=[]
for al in splitbam:
	chrom = splitbam.getrname(al.rname)
	read= al.qname
	start = al.pos
	end   = al.aend
	name  = splitbam.getrname(al.rnext)
	nstart = al.pnext
	rlen = al.rlen
	#print(chrom)
	#break
	if chrom=='%s'%options.chr:
		for key, value in myhier.items():
			if name in value: 
				name=key
				i = chrom, start, start+rlen, name, read, name
				myreads.append(i)
			elif chrom in value:
				chrom=key
				j = name, nstart, nstart+rlen, chrom, read, name
				myreads.append(j)
			else:
				pass

del(splitbam)
myreads=np.unique(myreads, axis=0)
myreads_pd=pd.DataFrame(myreads)
#print(myreads_pd.head(n=50))

tespd=pd.read_csv('%s/%s'%(options.dir, options.matrix), header=0)
tespd_col=tespd[["Chromosome","Start","End","TEfamily",'%s'%options.acc]]
tespd_col['Start']= tespd_col['Start'] - 1000
tespd_col[tespd_col['Start'] < 0] = 0
tespd_col['End']= tespd_col['End'] + 1000
tespd_col['TEfamily']=tespd_col['TEfamily'].apply(lambda x: str(x).replace('[','').replace(']',''))
tespd_col['TEfamily']=tespd_col['TEfamily'].apply(lambda x: str(x).replace("\'","")) 

myreads_bed = pybedtools.BedTool.from_dataframe(myreads_pd)
tes_bed = pybedtools.BedTool.from_dataframe(tespd_col)
#print(myreads_bed.head(n=50))
#print(tes_bed.head(n=50))
f1 = tes_bed.intersect(myreads_bed, wa=True, wb=True)
n1 = tes_bed.intersect(myreads_bed, v=True) 



#print(f1.head(n=50))
#print(n1.head(n=50))
#break

myoverlap=[]
for j in f1:
	counter=0
	chr1=j[0]
	posS=j[1]
	posE=j[2]
	te1=j[3]
	pa=j[4]
	chr2=j[5]
	rs=j[6]
	re=j[7]
	te2=j[8]
	mm='mismatch'
	result=chr1, posS, posE, te1, pa 
	mismatch = chr1, posS, posE, te1, mm
	if te1==te2 and chr1=='%s'%options.chr:
		myoverlap.append(result)
	elif te1!=te2 and chr1=='%s'%options.chr:
		myoverlap.append(mismatch)
	else:
		pass

myoverlap, my_count=np.unique(myoverlap,axis=0, return_counts=True)
myoverlap_pd=pd.DataFrame(myoverlap)
mycount_pd=pd.DataFrame(my_count)
#print(myoverlap_pd.head(n=20))
#print(mycount_pd.head(n=20))
#print(len(myoverlap_pd))
#print(len(mycount_pd))
frames = [myoverlap_pd, mycount_pd]
myresult = pd.concat(frames,ignore_index=True,axis=1)

#myresult.to_csv('ohman', sep='\t')
#print(myresult)
#break
#print(len(myresult))
myresult = myresult.drop_duplicates([0,1,2,3],keep= 'first')
del frames
del myoverlap_pd
del mycount_pd
del myoverlap
del my_count
#print(myresult)
#print(len(myresult))
#print(myresult.head(n=50))
#print(myresult.tail(50))


myproperpairs=[]
for read in pairbam:
	chrom = pairbam.getrname(read.rname)
	myname = read.qname
	start = read.pos
	end = read.aend
	if chrom=='%s'%options.chr:
		k= chrom, start , end, myname
		myproperpairs.append(k)
	else:
		pass
myproperpairs=np.unique(myproperpairs, axis=0)
mypairs_pd=pd.DataFrame(myproperpairs)
del(myproperpairs)
mypairs_bed = pybedtools.BedTool.from_dataframe(mypairs_pd)
del(mypairs_pd)
result_bed = pybedtools.BedTool.from_dataframe(myresult)
del(myresult)
#print(mypairs_bed.head(n=5))

f2=result_bed.intersect(mypairs_bed, wa=True, c=True)
#print(f2.head(n=50)) 

f3=n1.intersect(mypairs_bed, wa=True, c=True)
#print(f3.head(n=50))

out=open('%s.covfilt%s.bed'%(options.acc, options.chr), 'w')
for line in f2:
	if not line[0].startswith('C'):
		pass
	else:
		chr=line[0]
		start=int(line[1])+1000
		end=int(line[2])-1000
		te=line[3]
		pa=line[4]
		split=line[5]
		pair=line[6]
		if chr=='%s'%options.chr and pa!="mismatch" and int(split)>=3:
			out.write('%s,%s,%s,%s,1\n'%(chr,start,end,te))
		elif chr=='%s'%options.chr and pa!="mismatch" and int(split)<3 and int(pair)>=1:
			out.write('%s,%s,%s,%s,0\n'%(chr,start,end,te))
		elif chr=='%s'%options.chr and pa!="mismatch" and int(split)<3 and int(pair)<1:
			out.write('%s,%s,%s,%s,NA\n'%(chr,start,end,te))
		elif chr=='%s'%options.chr and pa=="mismatch" and int(pair)>=1:
			out.write('%s,%s,%s,%s,0\n'%(chr,start,end,te))
		elif chr=='%s'%options.chr and pa=="mismatch" and int(pair)<1:
			out.write('%s,%s,%s,%s,NA\n'%(chr,start,end,te))
		elif chr=='%s'%options.chr:
			out.write('%s,%s,%s,%s,NA\n'%(chr,start,end,te))
		else:
			pass

del f2
for line in f3:
	if not line[0].startswith('C'):
		pass
	else:
		chr=line[0]
		start=int(line[1])+1000
		end=int(line[2])-1000
		te=line[3]
		pa=line[4]
		pair=line[5]
		if chr=='%s'%options.chr and int(pair)>=1:
			out.write('%s,%s,%s,%s,0\n'%(chr,start,end,te))
		elif chr=='%s'%options.chr and int(pair)<1:
			out.write('%s,%s,%s,%s,NA\n'%(chr,start,end,te))
		elif chr=='%s'%options.chr:
			out.write('%s,%s,%s,%s,NA\n'%(chr,start,end,te))
		else:
			pass
out.close()
