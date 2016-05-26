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

import json,os,random,signal,sys

# Make pipeline friendly output for head and tail commands
signal.signal(signal.SIGPIPE,signal.SIG_DFL)

def tag2compliment(tag):
  conversions={'A':'T','a':'t','C':'G','c':'g','G':'C','g':'c','T':'A','t':'a','N':'N','n':'n'}
  return ''.join([conversions[i] for i in tag[::-1]])

if len(sys.argv)==1:
  print '''A collection of scripts to process fastq files. Usage:

   fastqstream.py breakdown [value] input.fastq
    Breakdown of [bases, lengths].

  fastqstream.py cleanstagger [stagger1,stagger2,...,staggerN] input.fastq > output.fastq
    Removes reads without stagger at begining and also clears away stagger.

  fastqstream.py has [string] input.fastq > match.fastq 2> nomatch.fastq
    Keep reads containing [string].

  fastqstream.py illcheck [string] input.fastq > output.fastq
    Keep reads containing [string] as the barcode identifier.

  fastqstream.py katrep [prefix] input.fastq > output.fastq
    Keep and trim reads with exact prefix.

  fastqstream.py katres [string] input.fastq > output.fastq
    Keep and trim reads with exact sequence.

  fastqstream.py lacks [string] input.fastq > nomatch.fastq 2> match.fastq
    Remove reads lacking [string].

  fastqstream.py maxlen [# nt] input.fastq > output.fastq
    Remove reads that are above a certain size.

  fastqstream.py minlen [# nt] input.fastq > output.fastq
    Remove reads that are below a certain size.

  fastqstream.py modhead [string] input.fastq > output.fastq
    Add additional text to headers.

  fastqstream.py qualfilt [phred] input.fastq > output.fastq
    Removes reads with overall mean quality less than [phred].

  fastqstream.py removeN input.fastq > output.fastq
    Remove N bases.

  fastqstream.py reverse input.fastq > output.fastq
    Reverse reads.

  fastqstream.py rmdup input.fastq > output.fastq
    Remove reads with duplicate headers.

  fastqstream.py rmheaders [header1,header2,...,headerN] input.fastq > output.fastq
    Remove specific reads.

  fastqstream.py split [string] input.fastq > left.fastq 2> right.fastq
    Split sequence into two pieces by [string].

  fastqstream.py strip [string] input.fastq > output.fastq
    Strip [string] and subsequent characters from sequences (end strip).

  fastqstream.py subset [%] input.fastq > output.fastq
    Randomly generate a subset of reads.

  fastqstream.py subset2 [%] input.fastq > output.fastq 2> remaining.fastq
    Randomly generate a subset of reads, allows for replacement.

  fastqstream.py tofasta input.fastq > output.fasta
    Convert to fasta format.

  fastqstream.py trimlen [# nt] input.fastq > output.fastq
    Shorten reads.
'''
  sys.exit(0)

command=sys.argv[1]
value=sys.argv[-2]

if sys.argv[-1]=='-':fid=sys.stdin
else:fid=open(sys.argv[-1])

if (command=='has')|(command=='lacks'):compliment=tag2compliment(value)

if command=='rmheaders':headers=set(value.split(','))
else:headers=set()

breakdown,longest={},0
while True:
  header=fid.readline().strip()
  if not len(header):break
  sequence=fid.readline().strip()
  x=fid.readline().strip()
  phred=fid.readline().strip()

  if command=='breakdown':
    if value=='lengths':
      i=len(sequence)
      if i not in breakdown:breakdown[i]=0
      breakdown[i]+=1

    if value=='bases':
      for n in sequence:
        if n not in breakdown:breakdown[n]=0
        breakdown[n]+=1

  if command=="cleanstagger":
    staggerset=sys.argv[2].upper().split(',')
    lss=len(staggerset)
    checkstagger=True
    stagcount=0
    while checkstagger==True:
      stagger=staggerset[stagcount]
      ls=len(stagger)
      stagcount+=1
      if stagcount==lss:checkstagger=False
      if sequence[:ls]==stagger:
        sys.stdout.write("%s\n%s\n%s\n%s\n"%(header,sequence[ls:],x,phred[ls:]))
        checkstagger=False

  if command=='has':
    if (sequence.find(value)>-1):sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))
    else:
      if (sequence.find(compliment)>-1):sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))
      else:sys.stderr.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))

  if command=="illcheck":
    barcode=header.split(':')[-1]
    if barcode==value:
      sys.stdout.write("%s\n%s\n%s\n%s\n"%(header,sequence,x,phred))

  if command=='katrep':
    pos=len(value)
    if value==sequence[:pos]:
      sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence[pos:],x,phred[pos:]))

  if command=="katres":
    p1=sequence.find(sys.argv[2])
    p2=len(sys.argv[2])
    if p1>=0:
      sys.stdout.write("%s\n%s\n%s\n%s\n"%(header,sequence[p1+p2:],x,phred[p1+p2:]))

  if command=='lacks':
    if (sequence.find(value)>-1):sys.stderr.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))
    else:
      if (sequence.find(compliment)>-1):sys.stderr.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))
      else:sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))

  if command=='maxlen':
    if len(sequence)<=int(value):
      sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))

  if command=='minlen':
    if len(sequence)>=int(value):
      sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))

  if command=='modhead':
    sys.stdout.write('%s\n%s\n%s\n%s\n' %(header+value,sequence,x,phred))

  if command=="qualfilt":
    value=sum([ord(i)-33 for i in phred])/float(len(phred))
    if value>int(sys.argv[2]):
      sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))

  if command=='removeN':
    ind=[i for i in xrange(len(sequence)) if sequence[i]!='N']
    sequence2=''.join([sequence[i] for i in ind])
    phred2=''.join([phred[i] for i in ind])
    sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence2,x,phred2))

  if command=='reverse':
      sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,tag2compliment(sequence),x,phred[::-1]))
  
  if command=='rmdup':
    if header not in headers:sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))
    else:headers.add(header)

  if command=='rmheaders':
    if header in headers:sys.stderr.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))
    else:sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))

  if command=='split':
    i=sequence.find(sys.argv[2])
    j=i+len(sys.argv[2])
    if (i>-1):
      sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence[:i],x,phred[:i]))
      sys.stderr.write('%s\n%s\n%s\n%s\n' %(header,sequence[j:],x,phred[j:]))

  if command=="strip":
    i=sequence.find(sys.argv[2])
    if i>=0:
      sequence,phred=sequence[:i],phred[:i]
    else:
      s=''.join([i for i in sys.argv[2][:-1]])
      loopCond=True
      while loopCond:
        if len(s)>3:
          if sequence[-len(s):]==s:
            sequence=sequence[:-len(s)]
            phred=phred[:-len(s)]
            loopCond=False
          else:s=s[:-1]
        else:loopCond=False

    sys.stdout.write("%s\n%s\n%s\n%s\n"%(header,sequence,x,phred))

  if command=='subset':
    r=random.random()
    if float(value)>r*100:
      sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))

  if command=='subset2':
    r=random.random()
    if float(value)>r*100:sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))
    else:sys.stderr.write('%s\n%s\n%s\n%s\n' %(header,sequence,x,phred))

  if command=='tofasta':sys.stdout.write('>%s\n%s\n' %(header[1:],sequence))

  if command=='trimlen':
    value=int(value)
    sys.stdout.write('%s\n%s\n%s\n%s\n' %(header,sequence[:value],x,phred[:value]))





if command=='breakdown':
  if value=='lengths':
    sys.stdout.write('#Read_length\tCount\n')
    for n in sorted(breakdown):
      sys.stdout.write('%s\t%s\n' %(n,breakdown[n]))

  if value=='bases':
    breakdown2={}
    for n in breakdown:
      m=n.upper()
      if m not in breakdown2:breakdown2[m]=0
      breakdown2[m]+=breakdown[n]

    sys.stdout.write('#Nucleotide\tCount\n')
    for n in sorted(breakdown2):
      sys.stdout.write('%s\t%s\n' %(n,breakdown2[n]))
