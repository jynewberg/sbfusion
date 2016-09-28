#!/usr/bin/env python

# Copyright (C) 2016 Justin Y. Newberg
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os,sys

# Standalone copies of library functions
## Get the reverse compliment of a sequence
def tag2compliment(tag):
  rosetta={'A':'T','a':'t','C':'G','c':'g','G':'C','g':'c','T':'A','t':'a','N':'N','n':'n'}
  return ''.join([rosetta[i] for i in tag[::-1]])


## Reduce sam flags to principal values
def sam_flag(x):
  binstr=bin(x)[2:][::-1]
  values=[2**i for i in xrange(len(binstr)) if binstr[i]=='1']
  return values



if len(sys.argv)==1:
  doc='''cigar usage:
  cigar M input.sam > aligned.fastq
  cigar S input.sam > unaligned.fastq
'''
  sys.stdout.write(doc)
  sys.exit(0)

command = sys.argv[1]

if sys.argv[-1]=='-':fid = sys.stdin
else:
  if os.path.exists(sys.argv[-1])==False:
    sys.stderr.write("Check input file\n")
    sys.exit(1)
  fid = open(sys.argv[-1],"rbU")

while True:
  line=fid.readline().strip('\n')
  if line=='':break
  if line[0]=='@':continue

  row=line.split('\t')
  header=row[0]
  sequence=row[9]
  _='+'
  phred=row[10]
  cigar=row[5]

  if cigar=='*':continue

  flags=sam_flag(int(row[1]))

  features,sequences,phreds,value,position = [],[],[],'',0
  st,en=0,0
  for i in cigar:
    o = ord(i)
    if (o>47) & (o<58): value+=i
    else:
      if i!='D':en=st+int(value)
      features.append([int(value),i])
      sequences.append(sequence[st:en])
      phreds.append(phred[st:en])
      value=''
      st=en

  i=0
  while i<len(features):
    if features[i][1]!='D':i+=1
    else:
      phreds.pop(i)
      sequences.pop(i)
      features.pop(i)

  i=0
  while i<len(features):
    if features[i][1]!='I':i+=1
    else:
      phreds[i-1]+=phreds.pop(i)
      sequences[i-1]+=sequences.pop(i)
      features[i-1][0]+=features.pop(i)[0]

  i=0
  while i<len(features)-1:
    if features[i][1]!=features[i+1][1]:i+=1
    else:
      sequences[i]+=sequences.pop(i+1)
      phreds[i]+=phreds.pop(i+1)
      features[i][0]+=features.pop(i+1)[0]

  maxes={'M':-1,'S':-1,'H':-1}
  argmaxes={'M':-1,'S':-1,'H':-1}
  index={'M':'','S':'','H':''}
  index2={'M':'','S':'','H':''}
  for i in xrange(len(features)):
    if len(index[features[i][1]])<len(sequences[i]):
      index[features[i][1]]=sequences[i]
      index2[features[i][1]]=phreds[i]
      maxes[features[i][1]]=len(sequences[i])
      argmaxes[features[i][1]]=i

  if command=='M':
    if len(index['M'])<20:continue
    read_side="left"
    if (argmaxes['S']<argmaxes['M'])&(maxes['S']>19):read_side="right"
    sequence=index['M']
    phred=index2['M']
    if 16 in flags:
      sequence=tag2compliment(sequence)
      phred=phred[::-1]
    sys.stdout.write("@%s\n%s\n%s\n%s\n" %(header,sequeunce,_,phred))
    sys.stderr.write("@%s\t%s\n" %(header,read_side))

  if command=='S':
    if len(index['S'])<20:continue
    read_side="left"
    if (argmaxes['M']<argmaxes['S'])&(maxes['M']>19):read_side="right"
    sequence=index['S']
    phred=index2['S']
    if 16 in flags:
      sequence=tag2compliment(sequence)
      phred=phred[::-1]
    sys.stdout.write("@%s\n%s\n%s\n%s\n" %(header,sequence,_,phred))
    sys.stderr.write("@%s\t%s\n" %(header,read_side))
