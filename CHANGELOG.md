# Update log

##2016-09-27
- Uploaded new version of cigar.py. The old version of cigar.py did account for the fact that bowtie2 reports reverse-compliment sequences for reads that align to negative strand. This is important to consider when a fusion read maps to a locus with two overlapping gene annotations that are anti-sense to each other.
- Example files have been uploaded accordingly.
