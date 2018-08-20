# exon-k-mer-extractor
Using sliding window to extract k-mers for annotation (specfically exons) in GTF format

## Dependencies:
1. Linux
2. Python 2.6

## Get Started:
### Input
python kmer_analysis.py -g {input_annotation_gtf}.gtf -k n -r m -o {output_bed}.bed

### Output BED Format
column 1: chromosome number (chr1, chr2, chr3...)
column 2: start location (1445550...)
column 3: end location (1445800...)
column 4: unique id (chr1-14403:14404:ED0...)

In column 4, there are four designations for relative position (IU, EU, ED, ID) followed by a number. 
IU = intron upstream
EU = exon upstream
followed by a number indicating how many nucleotides away from the start location

ED = exon downstream
ID = intron downstream
followed by a number indicating how many nucleotides away from the end location

