# README
Last modified 2016-09-27

##Introduction
The code accompanying this readme functions as a toolkit for identifying novel Sleeping Beauty fusion transcripts in RNA-seq data. 

##Dependencies and code installation
This code has been tested on Mac OS 10.10+ systems and RHEL 6.5+

This code requires Python 2.6.* or 2.7.* with numpy 1.4.* or 1.8.*, respectively. It also requires samtools 0.1.19-44428cd, bowtie2 2.2.5, and bedtools 2.17.0. 

While Python 2.* is simple to install manually, numpy can be quite a challenge to build properly (at least on a Mac). Since most Linux/Mac systems come with these two components pre-installed, it is assumed that you are using the system version of Python 2.7.* and numpy.

In order to build samtools, bowtie, and bedtools your system will need to have gcc (Linux) or XCode (Mac) installed. 

This toolbox contains an installation script (```install.sh```) for the remaining dependencies (samtools, bowtie, and bedtools). This script will create a directory, ```nnlab/sbfusion``` in your home directory that contains these pieces of software as well as custom python code and reference annotations.

Note that the reference data takes a few hours to index.

##The sbfusion workflow
The SB fusion workflow consists of two parts:

- identification of mouse-transposon fusions (from RNA-seq data), including mapping these fusion transcripts to corresponding sites in the mouse genome.
- utilization of these sites with insertions detected by SBCaptureSeq to define fusion transcripts. Then merging these novel transcripts into GTF files that can be used in cufflinks.

##Trying out the code
You can download a toy dataset from [Figshare](https://figshare.com/s/6338a2d4211140fe0c80). Make sure to unzip the file once you've gotten it (so the filename ends in fastq). 

###Running the first step (fusion detection)
You can run the first part of the sbfusion workflow as follows:

```
export PATH=~/nnlab/sbfusion/bin:$PATH

fastq2fusion.sh ~/nnlab/sbfusion/srv/mm9/vector.fa ~/nnlab/sbfusion/srv/mm9/onc2.fa sbfusion-toy.fastq
```

This takes about one minute to run on the toy dataset. 

###Running the second step (fusion transcript definition)
The second part of the workflow requires a transcript annotation file (in GTF format). The one we provide was derived from the refGene.txt.gz file from the [UCSC Genome Browser](http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database), downloaded on 2015-05-03. We converted the genePred file to GTF format using command line tools provided at that [site](http://hgdownload.cse.ucsc.edu/downloads.html#source_downloads). 

Additionally, this second step requires insertions defined by SBCaptureSeq. We have provided an example file that comes from the same tumor as the toy RNA-seq dataset (the installer places a copy of this bedfile in the ~/nnlab/sbfusion/var/example directory).

```
export PATH=~/nnlab/sbfusion/bin:$PATH

defineNewTranscripts.py ~/nnlab/sbfusion/srv/mm9/refGene.gtf ~/nnlab/sbfusion/var/example/exampleSBCapSeqInsertions.bed sbfusion-toy/sbfusion-toy.fusion.bed > fusionTranscripts.gtf
```

This takes about a minute to run.

This will produce a file, ```fusionTranscripts.gtf```, in your current working directory that contains transcript annotations in the gtf format. 

Note that in the ```data``` folder that accompanies this readme there are two files beginning in "REFERENCE_" that you can view to make sure your outputs match what is expected.
