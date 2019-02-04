from pathlib import Path
import subprocess
from shlex import split
import sys
"""
Needs samtools, bowtie2, cutadapt, fastqc in the PATH and this is in python3.7
"""

# Set up the dir to raw reads
path = "/Users/yzhang250/Desktop/ND-research/Characterization-of-tumor-derived-Microvesicles-TMV-microRNAs-across" \
	   "-different-cancer-cell-lines-s/Raw_reads/"

# This can be any SCDS fastq raw read files
fastq_file = ""
sample_name = ""
trimmed_fastq_file = ""

def extract_sample_name(fastq_file):
    return fastq_file[4:fastq_file.index("-")]

# Align index is default hsa_miRNA_mature_index, but can be hsa_miRNA_hairpin_index, hg19/hg19
align_index = "hsa_miRNA_mature_index"

# This step trimmed the NEB adapters
def cutadaptor():
	report = subprocess.check_output(['cutadapt',
							 '-f',
							 'fastq',
							 '-a',
							 'AGATCGGAAGAGCACACGTCT',
							 '-g',
							 'GTTCAGAGTTCTACAGTCCGACGATC',
							 '-e',
							 '0.15',
							 '-O',
							 '10',
							 '-m',
							 '14',
							 path + fastq_file,
							 '-o',
							 path + "trimmed_" + fastq_file])
	print("")
	print("============================================================================")
	print("Sample in analyzing: " + sample_name)
	print(report.decode("utf-8"))

# Make a miRNA hairpin and mature index file with Bowtie2
def make_mature_hairpin_bowtie2_index():
	hsa_miRNA_hairpin_index = Path(path + "hsa_miRNA_hairpin.fa")
	if not hsa_miRNA_hairpin_index.is_file():
		res = []
		hsa = False
		with open(path + "hairpin.fa") as file:
			for line in file:
				if line.startswith(">hsa"):
					hsa = True
					res.append(line)
				elif line.startswith(">"):
					hsa = False
				else:
					if hsa:
						res.append(line.replace("U", "T"))
		with open(path + "hsa_miRNA_hairpin.fa", "w") as file:
			for line in res:
				file.write(line)
		subprocess.check_output(["bowtie2-build", path + "hsa_miRNA_hairpin.fa", path + "hsa_miRNA_hairpin_index"])

	hsa_miRNA_mature_index = Path(path + "hsa_miRNA_mature.fa")
	if not hsa_miRNA_mature_index.is_file():
		res = []
		hsa = False
		with open(path + "mature.fa") as file:
			for line in file:
				if line.startswith(">hsa"):
					hsa = True
					res.append(line)
				elif line.startswith(">"):
					hsa = False
				else:
					if hsa:
						res.append(line.replace("U", "T"))
		with open(path + "hsa_miRNA_mature.fa", "w") as file:
			for line in res:
				file.write(line)
		subprocess.check_output(["bowtie2-build", path + "hsa_miRNA_mature.fa", path + "hsa_miRNA_mature_index"])

# Used bowtie2 to align the trimmed reads to align_index
# Method adapted from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4652620/
def align_to_index():
	subprocess.check_output(['bowtie2',
							 '-q',
							 '--phred33',
							 '-p',
							 '8',
							 '-D',
							 '20',
							 '--local',
							 '-R',
							 "3",
							 '-N',
							 '0',
							 '-L',
							 '8',
							 '-i',
							 "S,1,0.50",
							 '-x',
							 path + align_index,
							  path + trimmed_fastq_file,
							 '-S', path + "aln_" + sample_name +".sam"])
	print("")
	print("============================================================================" )
	print("Sample in analyzing: " + sample_name)
	print("Alignment is done! The result is in " + path + "aln_" + sample_name +".sam")

# Check the flatstat of the above aligned sam file
def check_sam_file():
	print("")
	print("============================================================================")
	print("Sample in analyzing: " + sample_name)
	report = subprocess.check_output(["samtools", "flagstat",  path + "aln_" + sample_name +".sam"])
	with open("Samtools_flagstat/" + sample_name + '_report.txt', "w") as outfile:
		for line in report.decode("utf-8"):
			outfile.write(line)


# Convert the sam file to bam file for downstream analysis
# Sort and index the bam file to make bai file
def convert_sort_index_sam_to_bam():
	my_cmd = ["samtools", "view", "-S", "-b", path + "aln_" + sample_name +".sam"]
	with open(path + "aln_" + sample_name +".bam", "w") as outfile:
		subprocess.call(my_cmd, stdout=outfile)
	subprocess.check_output(["samtools", "sort", "-o",
							 path + "sorted_" + "aln_" + sample_name +".bam", path + "aln_" + sample_name +".bam"])
	subprocess.check_output(["samtools", "index",
							 path + "sorted_" + "aln_" + sample_name +".bam"])

# Redirect the samtools idxstats to a file in Raw_counts
def generate_raw_counts():
	p1 = subprocess.Popen(split("samtools idxstats " + path + "sorted_" + "aln_" + sample_name + ".bam"),
						  stdout = subprocess.PIPE)
	with open("Raw_counts/" + sample_name + "_raw_count.txt", "w") as outfile:
		subprocess.Popen(split("cut -f1,3"), stdin = p1.stdout, stdout = outfile)

	print("")
	print("============================================================================")
	print("Sample in analyzing: " + sample_name)
	print("Raw counts file has been generated")

if __name__ == '__main__':
	fastq_files = sys.argv[1:]
	for file in fastq_files:
		fastq_file = file
		sample_name = extract_sample_name(fastq_file)
		trimmed_fastq_file = "trimmed_" + fastq_file
		cutadaptor()
		make_mature_hairpin_bowtie2_index()
		align_to_index()
		check_sam_file()
		convert_sort_index_sam_to_bam()
		generate_raw_counts()

