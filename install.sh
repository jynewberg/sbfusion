#!/usr/bin/env bash

initpath=`pwd`

rootdir="$HOME/nnlab/sbfusion"
execdir="$rootdir/bin"
srcdir="$rootdir/var/src"
buildir="$rootdir/var/build"
logdir="$rootdir/var/log"
srvdir="$rootdir/srv"
indexdir="$srvdir/mm9"

mkdir -p "$execdir" "$srcdir" "$buildir" "$logdir" "$indexdir"

export PATH=$execdir:$PATH

########################################################

toolbox="samtools"
version="0.1.19"
release="$toolbox"-"$version"
if [ ! -f "$srcdir"/"$release".tar.bz2 ]; then curl -L http://sourceforge.net/projects/samtools/files/"$toolbox"/"$version"/"$release".tar.bz2 > "$srcdir"/"$release".tar.bz2; fi
if [ ! -f "$execdir"/samtools ]; then
  rm -rf "$buildir"/"$release"
  tar -xf "$srcdir"/"$release".tar.bz2 -C "$buildir"
  cd "$buildir"/"$release"
  make > "$logdir"/"$toolbox".make.stdout 2> "$logdir"/"$toolbox".make.stderr
  cp samtools "$execdir"
fi

toolbox="bowtie2"
version="2.2.5"
release="$toolbox"-"$version"
if [ ! -f "$srcdir"/"$release"-source.zip ]; then curl -L http://sourceforge.net/projects/bowtie-bio/files/"$toolbox"/"$version"/"$release"-source.zip > "$srcdir"/"$release"-source.zip; fi
if [ ! -f "$execdir"/bowtie2 ]; then
  rm -rf "$buildir"/"$release"
  unzip -qq "$srcdir"/"$release"-source.zip -d "$buildir"
  cd "$buildir"/"$release"
  make > "$logdir"/"$toolbox".make.stdout 2> "$logdir"/"$toolbox".make.stderr
  cp bowtie2 bowtie2-align-l bowtie2-align-s bowtie2-build bowtie2-build-l bowtie2-build-s bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-s "$execdir"
fi

# note that 2.17.0 was originally used, but that was hosted on Google Code which has shut down. The newest version of bedtools is linked-to below. Bedtools has been quite stable for a few years, and light testing indices it still performs the same in this pipeline as before.
toolbox="bedtools2"
version="2.25.0"
release="$toolbox"-"$version"
if [ ! -f "$srcdir"/"$release".tar.gz ]; then curl -L curl -L https://github.com/arq5x/bedtools2/archive/v2.25.0.tar.gz > "$srcdir"/"$release".tar.gz; fi
if [ ! -f "$execdir"/bedtools ]; then
  rm -rf "$buildir"/"$release"
  tar -xf "$srcdir"/"$release".tar.gz -C "$buildir"
  cd "$buildir"/"$release"
  make > "$logdir"/"$toolbox".make.stdout 2> "$logdir"/"$toolbox".make.stderr
  cp bin/* "$execdir" > "$logdir"/"$toolbox".install.stdout 2> "$logdir"/"$toolbox".install.stderr
fi


########################################################

if [ ! -f "$indexdir"/refGene.gtf ];
then
  gzip -dc "$initpath"/data/refGene.gtf.gz > "$indexdir"/refGene.gtf
fi

if [ ! -f "$indexdir"/chromFa.tar.gz ]; 
then
  echo "Downloading mm9 reference"
  curl http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz > "$indexdir"/chromFa.tar.gz
fi

pylibdir=`python -c "import sys;print(sys.path[-1])"`
export PYTHONPATH=$pylibdir:$PYTHONPATH
mkdir -p $pylibdir


echo "Defining transposon vector sequence"
echo ">onc2
ccattcgccattcaggctgcgcaactgttgggaagggcgatcggtgcggg
cctcttcgctattacgccagctggcgaaagggggatgtgctgcaaggcga
ttaagttgggtaacgccagggttttcccagtcacgacgttgtaaaacgac
ggccagtgagcgcgcgtaatacgactcactatagggcgaattggagctcg
gatccctatacagttgaagtcggaagtttacatacacttaagttggagtc
attaaaactcgtttttcaactactccacaaatttcttgttaacaaacaat
agttttggcaagtcagttaggacatctactttgtgcatgacacaagtcat
ttttccaacaattgtttacagacagattatttcacttataattcactgta
tcacaattccagtgggtcagaagtttacatacactaagttgactgtgcct
ttaaacagcttggaaaattccagaaaatgatgtcatggctttagaagctt
gatggccgctctagaactaggattgcagcacgaaacaggaagctgactcc
acatggtcacatgctcactgaagtgttgacttccctgacagctgtgcact
ttctaaaccggttttctcattcatttacagttcagccgatgatgaaattg
ccgcactggttgttagcaacgtagccggtatgtgaaagatggattcgcgg
gaatttagtggatcccccgggctgcaggaattcgatctgaagcctataga
gtacgagccatagataaaataaaagattttatttagtctccagaaaaagg
ggggaatgaaagaccccacctgtaggtttggcaagctagcttaagtaacg
ccattttgcaaggcatggaaaatacataactgagaatagagaagttcaga
tcaaggttaggaacagagagacagcagaatatgggccaaacaggatatct
gtggtaagcagttcctgccccggctcagggccaagaacagatggtcccca
gatgcggtcccgccctcagcagtttctagagaaccatcagatgtttccag
ggtgccccaaggacctgaaaatgaccctgtgccttatttgaactaaccaa
tcagttcgcttctcgcttctgttcgcgcgcttctgctccccgagctcaat
aaaagagcccacaacccctcactcggcgcgccagtcctccgatagactgc
gtcgcccatcaagcttgctactagcaccagaacgcccgcgaggatctctc
aggtaataaagagcgccaaggctggctgcaagcggagcctctgagagcct
ctgagggccagggctactgcacccttggtcctcaacgctggggtcttcag
aactagaatgctgggggtggggtggggattcggttccctattccatcgcg
cgttaagatacattgatgagtttggacaaaccacaactagaatgcagtga
aaaaaatgctttatttgtgaaatttgtgatgctattgctttatttgtaac
cattataagctgcaataaacaagttggccgctcctgtgccagactctggc
gccgctgctctgtcaggtacctgttggtctgaaactcagccttgagcctc
tggagctgctcagcagtgaaggctgtgcgaggccgcttgtcctctttgtt
agggttcttcttctttggttttcgggacctgggacctggttgtcatggag
gagaaagggcagaggttactggttgctggagtctagctacttatccacaa
cccacgcacccaagcttgaggttgcagatactgggggtgggggggggggg
atgacccgcccaaggccatacaagtgttgggcattgggggtggtgatata
aacttgaggctgggcatgtgcccactgaccagaaggaaagtggtgtgtgt
gtgtgaaaatgagatggattggcagatgtagctaaaaggcctatcacaaa
ctaggggatctagcttgtggaaggctactcgaaatgtttgacccaagtta
aacaatttaaaggcaatgctaccaaatactaattgagtgtatgtaaactt
ctgacccactgggaatgtgatgaaagaaataaaagctgaaatgaatcatt
ctctctactattattctgatatttcacattcttaaaataaagtggtgatc
ctaactgacctaagacagggaatttttactaggattaaatgtcaggaatt
gtgaaaaagtgagtttaaatgtatttggctaaggtgtatgtaaacttccg
acttcaactgtatagggatcctctagctagagtcgacctcgagggggggc
ccggtacccagcttttgttccctttagtgagggttaatttcgagcttggc
gtaatcatggtcatagctgtttcctgtgtgaaattgttatccgctcacaa
ttccacacaacatacgagccggaagcataaagtgtaaagcctggggtgcc
taatgagtgagctaactcacattaattgcgttgcgctcactgcccgcttt
ccagtcgggaaacctgtcgtgccagctgcattaatgaatcggccaacgcg
cggggagaggcggtttgcgtattgggcgctcttccgcttcctcgctcact
gactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactc
aaaggcggtaatacggttatccacagaatcaggggataacgcaggaaaga
acatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcg
ttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaa
tcgacgctcaagtcagaggtggcgaaacccgacaggactataaagatacc
aggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctg
ccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgct
ttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgct
ccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgcc
ttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatc
gccactggcagcagccactggtaacaggattagcagagcgaggtatgtag
gcggtgctacagagttcttgaagtggtggcctaactacggctacactaga
aggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaa
aagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtg
gtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaa
gaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaa
ctcacgttaagggattttggtcatgagattatcaaaaaggatcttcacct
agatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatat
gagtaaacttggtctgacagttaccaatgcttaatcagtgaggcacctat
ctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcg
tgtagataactacgatacgggagggcttaccatctggccccagtgctgca
atgataccgcgagacccacgctcaccggctccagatttatcagcaataaa
ccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccg
cctccatccagtctattaattgttgccgggaagctagagtaagtagttcg
ccagttaatagtttgcgcaacgttgttgccattgctacaggcatcgtggt
gtcacgctcgtcgtttggtatggcttcattcagctccggttcccaacgat
caaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctcc
ttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcact
catggttatggcagcactgcataattctcttactgtcatgccatccgtaa
gatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatag
tgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataatac
cgcgccacatagcagaactttaaaagtgctcatcattggaaaacgttctt
cggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatg
taacccactcgtgcacccaactgatcttcagcatcttttactttcaccag
cgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaa
taagggcgacacggaaatgttgaatactcatactcttcctttttcaatat
tattgaagcatttatcagggttattgtctcatgagcggatacatatttga
atgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaa
aagtgccacctgacgcgccctgtagcggcgcattaagcgcggcgggtgtg
gtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgc
tcctttcgctttcttcccttcctttctcgccacgttcgccggctttcccc
gtcaagctctaaatcgggggctccctttagggttccgatttagtgcttta
cggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgg
gccatcgccctgatagacggtttttcgccctttgacgttggagtccacgt
tctttaatagtggactcttgttccaaactggaacaacactcaaccctatc
tcggtctattcttttgatttataagggattttgccgatttcggcctattg
gttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaa
tattaacgcttacaattt" > "$indexdir"/vector.fa


########################################################


echo "Installing python scripts"

cp "$initpath"/python/cigar.py "$execdir"
chmod 744 "$execdir"/cigar.py

cp "$initpath"/python/defineNewTranscripts.py "$execdir"
chmod 744 "$execdir"/defineNewTranscripts.py

cp "$initpath"/python/fastqstream.py "$execdir"
chmod 744 "$execdir"/fastqstream.py

cp "$initpath"/python/sam2fastq.py "$execdir"
chmod 744 "$execdir"/sam2fastq.py

cp "$initpath"/python/sbf_fmtbed.py "$execdir"
chmod 744 "$execdir"/sbf_fmtbed.py

echo "Installing bash scripts"

cp "$initpath"/bash/fastq2fusion.sh "$execdir"
chmod 744 "$execdir"/fastq2fusion.sh



########################################################

if [ ! -f "$indexdir"/onc2.fa ];
then
  echo "Making onc2.fa reference"
  tar -xf "$indexdir"/chromFa.tar.gz -O > "$indexdir"/onc2.fa

  cat "$indexdir"/vector.fa >> "$indexdir"/onc2.fa
fi

########################################################

if [ ! -f "$indexdir"/vector.1.bt2 ];
then
  echo "Generating Bowtie2 indexes on vector.fa"
  bowtie2-build "$indexdir"/vector.fa "$indexdir"/vector
fi

if [ ! -f "$indexdir"/onc2.1.bt2 ];
then
  echo "Generating Bowtie2 indexes on onc2.fa (takes 2 hours)"
  bowtie2-build "$indexdir"/onc2.fa "$indexdir"/onc2
fi

########################################################

# echo "Add example SBCaptureSeq insertion set that can be used with toy RNA-seq dataset
mkdir -p "$rootdir"/var/example
cp "$initpath"/data/exampleSBCapSeqInsertions.bed "$rootdir"/var/example

########################################################
