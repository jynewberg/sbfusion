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

import signal,sys

signal.signal(signal.SIGPIPE,signal.SIG_DFL)

# standalone copies of library functions
def sam_flag(x):
  binstr=bin(x)[2:][::-1]
  return [2**i for i in xrange(len(binstr)) if binstr[i]=='1' ]

def tag2compliment(tag):
  translator={'A':'T','a':'t','C':'G','c':'g','G':'C','g':'c','T':'A','t':'a','N':'N','n':'n'}
  return ''.join([translator[i] for i in tag[::-1]])

# Display documentation
if len(sys.argv)==1:
  doc="""Usage: sam2fastq.py [sequences.sam]
  Executable script to extract single end FASTQ data from a SAM file. 
  Dumps text to stdout.
"""
  sys.stdout.write(doc)
  sys.exit(0)


if sys.argv[-1]=='-':fid=sys.stdin
else:fid=open(sys.argv[-1],"rbU")

while True:
  line=fid.readline().strip('\n')
  if line=='':break
  if line[0]=='@':continue
  row=line.split('\t')

  header=row[0]
  flag=int(row[1])
  sequence=row[9]
  phred=row[10]

  if 16 in sam_flag(flag):
    sequence=tag2compliment(sequence)
    phred=phred[::-1]

  block="@%s\n%s\n+\n%s\n"%(header,sequence,phred)
  sys.stdout.write(block)

fid.close()
