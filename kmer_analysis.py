#python kmer_analysis.py -g /ysm-gpfs/home/zc264/scratch60/ENTex/rnaseqlib-clip/rnaseqlib_custom_annotation/gencode_basic/commonshortest/SE.hg38.basic.se.gtf -k 1 -r 300 -o /ysm-gpfs/home/zc264/scratch60/ENTex/kmer-analysis/SE.gencode.basic.1mer.bed
#python kmer_analysis.py -g /ysm-gpfs/home/zc264/scratch60/ENTex/miso/const_exons/gencode.v28.basic.annotation.min_10.const_exons.gtf -k 1 -r 300 -o /ysm-gpfs/home/zc264/scratch60/ENTex/kmer-analysis/CON.gencode.basic.1mer.bed

import csv
import subprocess
import argparse
import os

parser = argparse.ArgumentParser(description="k-mer extraction v0.1 - extract k-mer's in the flanking regions of annotations in GTF format.")

parser.add_argument('-g', '--gtf_file', type=str,
	help='A required annotation file')

parser.add_argument('-k', '--kmer_size', type=int,
	help='The size of each k-mer')

parser.add_argument('-r', '--flanking_range', type=int,
	help='Range of k-mer extraction (in bp) from the annotations')

parser.add_argument('-o', '--output_file', type=str,
	help='A required annotation file in bed format')

args = parser.parse_args()

if (not args.gtf_file) or (not os.path.isfile(args.gtf_file)):
	parser.error("Invalid GTF file, please specify an annotation in GTF format with the -g or --gtf_file option")

if (not args.kmer_size) or (not isinstance(args.kmer_size, (int, long))):
	parser.error("Invalid k-mer size, please specify a k-mer size in integer with the -k or --kmer_size option")

if (not args.flanking_range) or (not isinstance(args.flanking_range, (int, long))):
	parser.error("Invalid flanking range of k-mer extraction, please specify a range in integer with the -r or --flanking_range option")

if (not args.output_file):
	parser.error("Invalid output file, please specify an output directory with output file name with the -o or --output_file option")

K = args.kmer_size

#/ysm-gpfs/home/zc264/scratch60/ENTex/rnaseqlib-clip/rnaseqlib_custom_annotation/gencode_basic/commonshortest/SE.hg38.basic.se.gtf
#/ysm-gpfs/home/zc264/scratch60/ENTex/miso/const_exons/gencode.v28.basic.annotation.min_10.const_exons.gtf
gtf = None

with open(args.gtf_file,'r') as file1:
	reader = csv.reader(file1,delimiter='\t')
	gtf = [r for r in reader]
	file1.close()

with open(args.output_file,"w") as file2:
	print len(gtf)
	for exon in gtf:
		
		chromosome = exon[0]

		#find the start and end position of the annotation
		annotation_start = int(exon[3])
		annotation_end = int(exon[4])

		#find the length of the annotation with start and end nucleotide inclusive
		annotation_length = annotation_end - annotation_start + 1

		#find the length of half of the annotation
		annotation_length_half = (annotation_length-1)/2

		#start and end of k-mer frames
		anno_flank_start = annotation_start - K + 1
		anno_flank_end = annotation_end #annotation_end - 1 when indexed by range

		if exon[6] == "+":

			for j in range(anno_flank_start, anno_flank_end):
				if annotation_start + annotation_length_half > j + K: #if the end of the k-mer is less than the middle of the anntoation
					index = "EU" + str(abs(annotation_start - (j + K - 1)))
					file2.write("%s\t%s\t%s\t%s-%s:%s:%s\n" % (chromosome, j-1, j + K - 1, chromosome, j-1, j + K - 1, index))
				if annotation_end - annotation_length_half < j: # if the start of the k-mer is more than the middle of the annotation
					index = "ED" + str(abs(j -  annotation_end) - 1)
					file2.write("%s\t%s\t%s\t%s-%s:%s:%s\n" % (chromosome, j-1, j + K - 1, chromosome, j-1, j + K - 1, index))
				if (annotation_start + annotation_length_half == j + K - 1) or (annotation_end - annotation_length_half == j):
					if index == "EU" + str(annotation_length_half-1):
						index = "ED" + str(annotation_length_half-1)
					else:
						index = "EU" + str(annotation_length_half-1)
					file2.write("%s\t%s\t%s\t%s-%s:%s:%s\n" % (chromosome, j-1, j + K - 1, chromosome, j-1, j + K - 1, index))

			#analyze the 300 k-mers around the exon
			for i in range(1,args.flanking_range+1):

				#calculate the start and end position of the up and down stream k-mers in zero-base for python
				up_start = int(exon[3]) - i - 1
				up_stop = int(exon[3]) - i + K - 1
				down_start = int(exon[4]) + i - K + 1
				down_stop = int(exon[4]) + i + 1

				file2.write("%s\t%s\t%s\t%s-%s:%s:%s\n" % (chromosome, up_start, up_stop, chromosome, up_start, up_stop, "IU"+str(i),))
				file2.write("%s\t%s\t%s\t%s-%s:%s:%s\n" % (chromosome, down_start, down_stop, chromosome, down_start, down_stop, "ID"+str(i),))

		if exon[6] == "-":

			for j in range(anno_flank_start, anno_flank_end):
				if annotation_start + annotation_length_half > j + K: #if the end of the k-mer is less than the middle of the anntoation
					index = "ED" + str(abs(annotation_start - (j + K - 1)))
					#subprocess.call("echo %s$'\t'%s$'\t'%s$'\t' %s-%s:%s:%s >> %s"
					#	% (chromosome, j, j + K, chromosome, j, j + K, index, args.output_file), shell=True)
					file2.write("%s\t%s\t%s\t%s-%s:%s:%s\n" % (chromosome, j-1, j + K - 1, chromosome, j-1, j + K - 1, index))
				if annotation_end - annotation_length_half < j: # if the start of the k-mer is more than the middle of the annotation
					index = "EU" + str(abs(j -  annotation_end) - 1)
					file2.write("%s\t%s\t%s\t%s-%s:%s:%s\n" % (chromosome, j-1, j + K - 1, chromosome, j-1, j + K - 1, index))
				if (annotation_start + annotation_length_half == j + K - 1) or (annotation_end - annotation_length_half == j):
					if index == "ED" + str(annotation_length_half-1):
						index = "EU" + str(annotation_length_half-1)
					else:
						index = "ED" + str(annotation_length_half-1)
					file2.write("%s\t%s\t%s\t%s-%s:%s:%s\n" % (chromosome, j-1, j + K - 1, chromosome, j-1, j + K - 1, index))
			
			#analyze the 300 k-mers around the exon
			for i in range(1,args.flanking_range+1):

				#calculate the start and end position of the up and down stream k-mers in zero-base for python
				up_start = int(exon[3]) - i - 1
				up_stop = int(exon[3]) - i + K - 1
				down_start = int(exon[4]) + i - K + 1
				down_stop = int(exon[4]) + i + 1

				file2.write("%s\t%s\t%s\t%s-%s:%s:%s\n" % (chromosome, up_start, up_stop, chromosome, up_start, up_stop, "ID"+str(i),))
				file2.write("%s\t%s\t%s\t%s-%s:%s:%s\n" % (chromosome, down_start, down_stop, chromosome, down_start, down_stop, "IU"+str(i),))
