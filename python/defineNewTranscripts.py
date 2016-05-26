#!/usr/bin/env python

import marshal,os,sys
from numpy import array

# allow padding
def photrans(exonStart,exonEnd,intronStart,intronEnd,location):
  index={}
  for i in xrange(len(exonStart)):
    exon=range(exonStart[i]-2,exonEnd[i]+1+2)
    for j in exon:index[j]="Ex%s" %(i)

  for i in xrange(len(intronStart)):
    intron=range(intronStart[i],intronEnd[i]+1)
    for j in intron:
      if j not in index:index[j]="In%s" %(i+1)

  return "In0" if location not in index else index[location]

def run_cmd(x,y=None,verbose=True):
  import glob,shlex,sys
  from subprocess import Popen,PIPE

  for cmd in x.split('|'):
    args=shlex.split(cmd)
    pipe=Popen(args,stdout=PIPE, stdin=PIPE, stderr=PIPE)
    y,err=pipe.communicate(input=y)
    if (len(err)>0) & verbose:sys.stderr.write(err+'\n')

  return y

def argsort(seq):
  return sorted(range(len(seq)),key=seq.__getitem__)

def find(statement):
  from numpy import nonzero, ravel
  return nonzero(ravel(statement))[0]











# Parse user inputs
if sys.argv[1]=='-':gtf_file=sys.stdin
else:gtf_file=open(sys.argv[1],"rbU")

# Set hard-coded parameters
offset=15000
offset_type="abs"

# Define transcripts from the GTF file
index_transcript={}
while True:
  line=gtf_file.readline().strip('\n')
  if line=='':break
  row=line.split('\t')
  fields={field.split(' ')[0]:field.split(' ')[1][1:-1] for field in row[8].strip(';').split("; ")}

  transcript=fields["transcript_id"]
  gene='' if "gene_name" not in fields else fields["gene_name"]

  if transcript not in index_transcript:index_transcript[transcript]={"chrom":row[0],"strand":row[6],"gene":gene,"transcript":transcript,"overlapping":set(),"locus":''}

  feature=row[2]
  if feature+"Start" not in index_transcript[transcript]:index_transcript[transcript][feature+"Start"]=[]
  if feature+"End" not in index_transcript[transcript]:index_transcript[transcript][feature+"End"]=[]
  index_transcript[transcript][feature+"Start"].append(int(row[3]))
  index_transcript[transcript][feature+"End"].append(int(row[4]))

# Sort exons and define introns
for transcript in index_transcript:
  ind=argsort([i for i in index_transcript[transcript]["exonStart"]])
  index_transcript[transcript]["exonStart"]=[index_transcript[transcript]["exonStart"][i] for i in ind]
  index_transcript[transcript]["exonEnd"]=[index_transcript[transcript]["exonEnd"][i] for i in ind]
  index_transcript[transcript]["txStart"]=index_transcript[transcript]["exonStart"][0]
  index_transcript[transcript]["txEnd"]=index_transcript[transcript]["exonEnd"][-1]

  inStart=index_transcript[transcript]["exonStart"][0]+1
  inEnd=-1

  if "intronStart" not in index_transcript[transcript]:index_transcript[transcript]["intronStart"]=[]
  if "intronEnd" not in index_transcript[transcript]:index_transcript[transcript]["intronEnd"]=[]

  for i in xrange(len(index_transcript[transcript]["exonStart"])-1):
    inEnd=index_transcript[transcript]["exonStart"][i+1]-1
    index_transcript[transcript]["intronStart"].append(inStart)
    index_transcript[transcript]["intronEnd"].append(inEnd)
    inStart=index_transcript[transcript]["exonEnd"][i+1]+1

  index_transcript[transcript]["locus"]="%s:%s-%s" %(index_transcript[transcript]["chrom"],index_transcript[transcript]["exonStart"][0],index_transcript[transcript]["exonEnd"][-1])

  if index_transcript[transcript]["strand"]=='-':
    if "CDSStart" in index_transcript[transcript]:
      index_transcript[transcript]["CDSStart"]=index_transcript[transcript]["CDSStart"][::-1]
      index_transcript[transcript]["CDSEnd"]=index_transcript[transcript]["CDSEnd"][::-1]
    index_transcript[transcript]["exonStart"]=index_transcript[transcript]["exonStart"][::-1]
    index_transcript[transcript]["exonEnd"]=index_transcript[transcript]["exonEnd"][::-1]
    index_transcript[transcript]["intronStart"]=index_transcript[transcript]["intronStart"][::-1]
    index_transcript[transcript]["intronEnd"]=index_transcript[transcript]["intronEnd"][::-1]

# Identify overlapping transcripts
index_gene={}
for transcript in index_transcript:
  gene=index_transcript[transcript]["gene"]
  if gene not in index_gene:index_gene[gene]=set()
  index_gene[gene].add(transcript)

for transcript in index_transcript:
  gene=index_transcript[transcript]["gene"]
  if gene!='':
    transcripts=index_gene[gene].copy()
    transcripts.remove(transcript)
    index_transcript[transcript]["overlapping"]=transcripts

if offset_type=='abs':offset_value=offset

# Index transcript start and stop sites on each chromosome to speed up annotation of BEDfile-like data. 
index_chrom={}
for i in index_transcript:
  chrom=index_transcript[i]["chrom"]
  txStart=index_transcript[i]["txStart"]
  txEnd=index_transcript[i]["txEnd"]
  strand=index_transcript[i]["strand"]

  if offset_type=='pct':offset_value=int((txEnd-txStart)*offset/100.+0.5)
  adjStart=txStart-offset_value*(strand=='+')
  adjEnd=txEnd+offset_value*(strand=='-')

  if chrom not in index_chrom:index_chrom[chrom]={"adjStart":[],"adjEnd":[],"txStart":[],"txEnd":[],"transcript":[]}
  index_chrom[chrom]["txStart"].append(txStart)
  index_chrom[chrom]["txEnd"].append(txEnd)
  index_chrom[chrom]["adjStart"].append(adjStart)
  index_chrom[chrom]["adjEnd"].append(adjEnd)
  index_chrom[chrom]["transcript"].append(index_transcript[i]["transcript"])












capfile=sys.argv[2] #"data/SBC_D5SP_R001.sbc.bed"
if sys.argv[3]=='-':fusfid=sys.stdin
else:fusfid=open(sys.argv[3],"rbU")

#index=marshal.load(open(marfile,"rbU"))
#index_transcript=index["transcript"]
#index_chrom=index["chrom"]

# Make the lists in index_chrom arrays for faster searching. Not done previously so the index_gtf script wouldn't have to load numpy. Probably should have been done then, though.
for i in index_chrom: 
  index_chrom[i]["adjStart"]=array(index_chrom[i]["adjStart"])
  index_chrom[i]["adjEnd"]=array(index_chrom[i]["adjEnd"])
  index_chrom[i]["txStart"]=array(index_chrom[i]["txStart"])
  index_chrom[i]["txEnd"]=array(index_chrom[i]["txEnd"])
  index_chrom[i]["transcript"]=array(index_chrom[i]["transcript"])
  
insertions=[line.split('\t') for line in open(capfile,"rbU").read().strip('\n').split('\n')]
for chrom,chromStart,chromEnd,name,score,strand in insertions:
  chromStart,chromEnd,score=int(chromStart),int(chromEnd),int(score)

  if chromEnd==chromStart+1:idx=find((chromStart>=index_chrom[chrom]["adjStart"]) & (chromStart<index_chrom[chrom]["adjEnd"])) if chrom in index_chrom else []
  else:idx=find(((chromStart>=index_chrom[chrom]["adjStart"]) & (chromStart<=index_chrom[chrom]["adjEnd"])) | ((index_chrom[chrom]["adjStart"]>=chromStart) & (index_chrom[chrom]["adjStart"]<=chromEnd)) | ((index_chrom[chrom]["adjEnd"]>=chromStart) & (index_chrom[chrom]["adjEnd"]<=chromEnd))) if chrom in index_chrom else []

  transcripts=sorted(set([index_chrom[chrom]["transcript"][i] for i in idx]))

  for transcript in transcripts:
    if "insertion" not in index_transcript[transcript]:index_transcript[transcript]["insertion"]={}
    if strand not in index_transcript[transcript]["insertion"]:index_transcript[transcript]["insertion"][strand]={}
    if chromStart not in index_transcript[transcript]["insertion"][strand]:index_transcript[transcript]["insertion"][strand][chromStart]=0
    index_transcript[transcript]["insertion"][strand][chromStart]+=score

fusions=[]
while True:
  line=fusfid.readline().strip('\n')
  if line=='':break
  fusions.append(line.split('\t'))

#fusions=[line.split('\t') for line in open("data/WTS_D5SP.fusion.bed","rbU").read().strip('\n').split('\n')]
#fusions+=[line.split('\t') for line in open("data/mRNA_D5SP_R001.fusion.bed","rbU").read().strip('\n').split('\n')]
for chrom,chromStart,chromEnd,name,score,strand in fusions:
  chromStart,chromEnd,score=int(chromStart),int(chromEnd),int(score)

  if chromEnd==chromStart+1:idx=find((chromStart>=index_chrom[chrom]["adjStart"]) & (chromStart<index_chrom[chrom]["adjEnd"])) if chrom in index_chrom else []
  else:idx=find(((chromStart>=index_chrom[chrom]["adjStart"]) & (chromStart<=index_chrom[chrom]["adjEnd"])) | ((index_chrom[chrom]["adjStart"]>=chromStart) & (index_chrom[chrom]["adjStart"]<=chromEnd)) | ((index_chrom[chrom]["adjEnd"]>=chromStart) & (index_chrom[chrom]["adjEnd"]<=chromEnd))) if chrom in index_chrom else []
  #if chromEnd==chromStart+1:idx=find((chromStart>=index_chrom[chrom]["txStart"]) & (chromStart<index_chrom[chrom]["txEnd"])) if chrom in index_chrom else []
  #else:idx=find(((chromStart>=index_chrom[chrom]["txStart"]) & (chromStart<=index_chrom[chrom]["txEnd"])) | ((index_chrom[chrom]["txStart"]>=chromStart) & (index_chrom[chrom]["txStart"]<=chromEnd)) | ((index_chrom[chrom]["txEnd"]>=chromStart) & (index_chrom[chrom]["txEnd"]<=chromEnd))) if chrom in index_chrom else []

  transcripts=sorted(set([index_chrom[chrom]["transcript"][i] for i in idx]))

  for transcript in transcripts:
    if "fusion" not in index_transcript[transcript]:index_transcript[transcript]["fusion"]={}
    if strand not in index_transcript[transcript]["fusion"]:index_transcript[transcript]["fusion"][strand]={}
    if name not in index_transcript[transcript]["fusion"][strand]:index_transcript[transcript]["fusion"][strand][name]={}
    for i in range(chromStart,chromEnd+1):
      if i not in index_transcript[transcript]["fusion"][strand][name]:index_transcript[transcript]["fusion"][strand][name][i]=0
      index_transcript[transcript]["fusion"][strand][name][i]+=score

# merge fusions into sequences
for transcript in index_transcript:
  if "fusion" not in index_transcript[transcript]:continue
  fusions=index_transcript[transcript]["fusion"]
  for strand in fusions:
    for name in fusions[strand]:
      merged={}
      i=sorted(fusions[strand][name])
      st,la=i[0],i[0]
      running_score=0
      for j in i[1:]:
        if j==la+1:
          if fusions[strand][name][j]>running_score:
            running_score=fusions[strand][name][j]
          la=j
        if j>la+1:
          key="%s-%s" %(st,la)
          merged[key]=running_score
          running_score=0
          st,la=j,j

      key="%s-%s" %(st,j)
      merged[key]=running_score
      fusions[strand][name]=merged

# annotate insertions
for transcript in index_transcript:
  if "insertion" not in index_transcript[transcript]:continue
  strand=index_transcript[transcript]["strand"]
  insertions=index_transcript[transcript]["insertion"]
  for orientation in insertions:
    if orientation==strand:
      if "activating_insertion" not in index_transcript[transcript]:index_transcript[transcript]["activating_insertion"]={}
      index_transcript[transcript]["activating_insertion"]=insertions[orientation]
    if orientation!=strand:
      if "inactivating_insertion" not in index_transcript[transcript]:index_transcript[transcript]["inactivating_insertion"]={}
      index_transcript[transcript]["inactivating_insertion"]=insertions[orientation]

# annotate fusions
for transcript in index_transcript:
  if "fusion" not in index_transcript[transcript]:continue
  strand=index_transcript[transcript]["strand"]
  antisense='+' if strand=='-' else '-'

  index_transcript[transcript]["cutsite"]={"LunSD":set(),"SA_exon":set(),"exon_En2SA":set()}

  exon_bases=[]
  for i in xrange(len(index_transcript[transcript]["exonStart"])):
    exon=range(index_transcript[transcript]["exonStart"][i],index_transcript[transcript]["exonEnd"][i]+1)
    exon_bases+=exon

  exon_bases=set(exon_bases)

  if "activating_insertion" in index_transcript[transcript]:
    insertion_sites=[insertion_site for insertion_site in index_transcript[transcript]["activating_insertion"]]
    insertion_site_down=min(insertion_sites) if strand=="+" else max(insertion_sites)
    insertion_site_up=max(insertion_sites) if strand=="-" else max(insertion_sites)

    fusions=index_transcript[transcript]["fusion"]

    if strand in fusions:
      for name in fusions[strand]:
        for fusion_site in fusions[strand][name]:
          fspl=fusion_site.split('-')
          fleft,fright=int(fspl[0]),int(fspl[1])

          fusion_bases=set(range(fleft,fright+1))
          if len(fusion_bases.intersection(exon_bases))==0:continue

          if name=="LunSD":
            if (insertion_site_down-fleft<0)&(strand=='+'):
              if "sb_splice_donor" not in index_transcript[transcript]:index_transcript[transcript]["sb_splice_donor"]={}
              index_transcript[transcript]["sb_splice_donor"][fusion_site]=fusions[strand]["LunSD"][fusion_site]
              pt=photrans(index_transcript[transcript]["exonStart"],index_transcript[transcript]["exonEnd"],index_transcript[transcript]["intronStart"],index_transcript[transcript]["intronEnd"],fleft)
              index_transcript[transcript]["cutsite"]["LunSD"].add(pt)
            if (insertion_site_down-fright>0)&(strand=='-'):
              if "sb_splice_donor" not in index_transcript[transcript]:index_transcript[transcript]["sb_splice_donor"]={}
              index_transcript[transcript]["sb_splice_donor"][fusion_site]=fusions[strand]["LunSD"][fusion_site]
              pt=photrans(index_transcript[transcript]["exonStart"],index_transcript[transcript]["exonEnd"],index_transcript[transcript]["intronStart"],index_transcript[transcript]["intronEnd"],fright)
              index_transcript[transcript]["cutsite"]["LunSD"].add(pt)

          if name=="SA_exon":
            if (fright-insertion_site_up<0)&(strand=='+'):
              if "sb_splice_acceptor" not in index_transcript[transcript]:index_transcript[transcript]["sb_splice_acceptor"]={}
              index_transcript[transcript]["sb_splice_acceptor"][fusion_site]=fusions[strand]["SA_exon"][fusion_site]
              pt=photrans(index_transcript[transcript]["exonStart"],index_transcript[transcript]["exonEnd"],index_transcript[transcript]["intronStart"],index_transcript[transcript]["intronEnd"],fright)
              index_transcript[transcript]["cutsite"]["SA_exon"].add(pt)
            if (fleft-insertion_site_up>0)&(strand=='-'):
              if "sb_splice_acceptor" not in index_transcript[transcript]:index_transcript[transcript]["sb_splice_acceptor"]={}
              index_transcript[transcript]["sb_splice_acceptor"][fusion_site]=fusions[strand]["SA_exon"][fusion_site]
              pt=photrans(index_transcript[transcript]["exonStart"],index_transcript[transcript]["exonEnd"],index_transcript[transcript]["intronStart"],index_transcript[transcript]["intronEnd"],fleft)
              index_transcript[transcript]["cutsite"]["SA_exon"].add(pt)

  if "inactivating_insertion" in index_transcript[transcript]:
    insertion_sites=[insertion_site for insertion_site in index_transcript[transcript]["inactivating_insertion"]]
    insertion_site_down=max(insertion_sites) if strand=="+" else max(insertion_sites)
    insertion_site_up=min(insertion_sites) if strand=="-" else max(insertion_sites)

    fusions=index_transcript[transcript]["fusion"]

    if antisense in fusions:
      for name in fusions[antisense]:
        for fusion_site in fusions[antisense][name]:
          fspl=fusion_site.split('-')
          fleft,fright=int(fspl[0]),int(fspl[1])

          fusion_bases=set(range(fleft,fright+1))
          if len(fusion_bases.intersection(exon_bases))==0:continue

          if name=="exon_En2SA":
            if (insertion_site_down-fright>0)&(strand=='+'):
              if "neo_splice_acceptor" not in index_transcript[transcript]:index_transcript[transcript]["neo_splice_acceptor"]={}
              index_transcript[transcript]["neo_splice_acceptor"][fusion_site]=fusions[antisense]["exon_En2SA"][fusion_site]
              pt=photrans(index_transcript[transcript]["exonStart"],index_transcript[transcript]["exonEnd"],index_transcript[transcript]["intronStart"],index_transcript[transcript]["intronEnd"],fright)
              index_transcript[transcript]["cutsite"]["exon_En2SA"].add(pt)
            if (insertion_site_down-fleft<0)&(strand=='-'):
              if "neo_splice_acceptor" not in index_transcript[transcript]:index_transcript[transcript]["neo_splice_acceptor"]={}
              index_transcript[transcript]["neo_splice_acceptor"][fusion_site]=fusions[antisense]["exon_En2SA"][fusion_site]
              pt=photrans(index_transcript[transcript]["exonStart"],index_transcript[transcript]["exonEnd"],index_transcript[transcript]["intronStart"],index_transcript[transcript]["intronEnd"],fleft)
              index_transcript[transcript]["cutsite"]["exon_En2SA"].add(pt)


for transcript in index_transcript:
  chrom=index_transcript[transcript]["chrom"]
  strand=index_transcript[transcript]["strand"]
  gene=index_transcript[transcript]["gene"]
  if "activating_insertion" in index_transcript[transcript]:
    if "sb_splice_acceptor" in index_transcript[transcript]:
      for cutsite in sorted(set([int(i[2:]) for i in index_transcript[transcript]["cutsite"]["SA_exon"]])):
        if cutsite+1>0:
          exons_upstream=range(0,cutsite+1)
          workex=0
          for i in exons_upstream:
            transcript2=transcript.replace("NM_","SB_exonSA_%s_" %(cutsite+1))
            workex+=1
            chromStart=index_transcript[transcript]["exonStart"][i]
            chromEnd=index_transcript[transcript]["exonEnd"][i]
            line="%s\tstdin\texon\t%s\t%s\t.\t%s\t.\t" %(chrom,chromStart,chromEnd,strand)
            line+="""gene_id "%s"; transcript_id "%s"; exon_number "%s"; exon_id "%s"; gene_name "%s";""" %(gene,transcript2,workex,transcript2+'.'+str(workex),gene)
            print line

for transcript in index_transcript:
  chrom=index_transcript[transcript]["chrom"]
  strand=index_transcript[transcript]["strand"]
  gene=index_transcript[transcript]["gene"]
  if "activating_insertion" in index_transcript[transcript]:
    if "sb_splice_donor" in index_transcript[transcript]:
      for cutsite in sorted(set([int(i[2:]) for i in index_transcript[transcript]["cutsite"]["LunSD"]])):
        if cutsite<len(index_transcript[transcript]["exonEnd"]):
          exons_downstream=range(cutsite,len(index_transcript[transcript]["exonEnd"]))
          workex=0
          for i in exons_downstream:
            transcript2=transcript.replace("NM_","SB_LunSD_%s_" %cutsite)
            workex+=1
            chromStart=index_transcript[transcript]["exonStart"][i]
            chromEnd=index_transcript[transcript]["exonEnd"][i]
            line="%s\tstdin\texon\t%s\t%s\t.\t%s\t.\t" %(chrom,chromStart,chromEnd,strand)
            line+="""gene_id "%s"; transcript_id "%s"; exon_number "%s"; exon_id "%s"; gene_name "%s";""" %(gene,transcript2,workex,transcript2+'.'+str(workex),gene)
            print line

for transcript in index_transcript:
  chrom=index_transcript[transcript]["chrom"]
  strand=index_transcript[transcript]["strand"]
  gene=index_transcript[transcript]["gene"]
  if "inactivating_insertion" in index_transcript[transcript]:
    if ("neo_splice_acceptor" in index_transcript[transcript])|("neo_splice_acceptor" in index_transcript[transcript]):
      for cutsite in sorted(set([int(i[2:]) for i in index_transcript[transcript]["cutsite"]["exon_En2SA"]])):
        if cutsite+1<len(index_transcript[transcript]["exonEnd"]):
          exons_upstream=range(0,cutsite+1)
          workex=0
          for i in exons_upstream:
            transcript2=transcript.replace("NM_","SB_En2SA_%s_" %(cutsite+1))
            workex+=1
            chromStart=index_transcript[transcript]["exonStart"][i]
            chromEnd=index_transcript[transcript]["exonEnd"][i]
            line="%s\tstdin\texon\t%s\t%s\t.\t%s\t.\t" %(chrom,chromStart,chromEnd,strand)
            line+="""gene_id "%s"; transcript_id "%s"; exon_number "%s"; exon_id "%s"; gene_name "%s";""" %(gene,transcript2,workex,transcript2+'.'+str(workex),gene)
            print line
