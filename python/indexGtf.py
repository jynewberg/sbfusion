#!/usr/bin/env python

import marshal,sys


# Purely pythonic way to sort lists. Requires no external toolboxes.
# Copied from library.py. 
def argsort(seq):
  return sorted(range(len(seq)),key=seq.__getitem__)


# Parse user inputs
if sys.argv[-1]=='-':gtf_file=sys.stdin
else:gtf_file=open(sys.argv[-1],"rbU")
output_path=sys.argv[-2]

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

# Dump indexes to a file
index={}
index['transcript']=index_transcript
index['chrom']=index_chrom
marshal.dump(index,open(output_path,"wb"))
