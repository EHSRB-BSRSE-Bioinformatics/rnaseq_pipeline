import re
import pandas as pd
from glob import glob
import pathlib

config_file_name = "config.yaml"

configfile: config_file_name

print("config file: " + config_file_name)

# set up conditions and samples
sample_info = config["samples"]

SAMPLES = sum([i.split(" ") for i in list(sample_info.values())], [])
CONDITIONS = list(sample_info.keys())

print("conditions: " + str(CONDITIONS))
print("samples: " + str(SAMPLES))

# set up directories
root = Path("/")

genome_dir = root / config["genome_root_dir"] / config["genome"]

input_dir = Path(config["input_dir"])
output_dir = Path(config["output_dir"])

trim_dir = output_dir / "trimmed"
fastqc_dir = output_dir / "fastqc"
align_dir = output_dir / "aligned"
procesed_dir = output_dir / "processed"
log_dir = output_dir / "logs"

output_dir.mkdir(parents=True, exist_ok=True)
trim_dir.mkdir(parents=True, exist_ok=True)
fastqc_dir.mkdir(parents=True, exist_ok=True)
align_dir.mkdir(parents=True, exist_ok=True)
procesed_dir.mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)

print("using genome: " + str(genome_dir))

print("input dir: " + str(input_dir))
print("output dir: " + str(output_dir))

print("trimmed sequences dir: " + str(trim_dir))
print("fastqc dir: " + str(fastqc_dir))
print("alignments dir: " + str(align_dir))
print("snakemake log dir: " + str(log_dir))

shell("cp -rf Snakefile {log_dir}")
shell("cp -rf {config_file_name} {log_dir}")

# gene names file, needed for annotation
# TODO do we need this?

rule gene_names:
    input: genome_dir / config["annotation_filename"]
    output: genome_dir / "geneid_to_name.txt"
    run:
        shell("zcat {input} | awk 'BEGIN{FS=\"\t\"}{split($9,a,\";\"); if($3~\"gene\") print) a[1]\"\t\"a[3]}' | sed 's/gene_source \"ensembl\"//' | tr ' ' '\t' | cut -f 2,5 | tr -d '\"' > {output}")

# make a metadata table from the sample names
rule metadata:
    input:
        expand(str(input_dir / "{sample}.fastq.gz"), sample=SAMPLES)
    output:
        output_dir / "metadata.csv"
    run:
        regex_use = config.regex_metadata.replace("{","{{",-1).replace("}","}}",-1) #this suppresses wildcard expansion by snakemake. Otherwise, snakemake thinks there are wildcards in the regular expression.
        if config.first_in_regex:
            shell( "python Lib/SampleMetadata.py -i %s/ -o %s -r '%s'" %(indir, output, regex_use) )
        else:
            shell("python Lib/SampleMetadata.py -i %s/ -o %s -r '%s' -f" %(indir,output, regex_use))

# run whole pipeline

rule all:
    input:
        expand(str(trim_dir / "{sample}_trimmed.fq.gz"), sample=SAMPLES),
        expand(str(align_dir / "{sample}_trimmed_bismark_bt2.bam"), sample=SAMPLES),
        output_dir / "multiqc_report.html",
    log: log_dir / "all.log"

##################################
### Trimming raw reads : Fastp ###
##################################

################################
#INFORMATION ON TRIMMING PROCESS
#--cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 == Trimmomatic  "LEADING:3"
#--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 == Trimmomatic  "TRAILING:3"
#--cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 == Trimmomatic  "SLIDINGWINDOW:4:15"

#### If additional trimming is needed (see multiQCreport): 
# add to R1: --trim_front1 {amount_bases} and/or --trim_tail1 {amount_bases}
# add to R2: --trim_front2 {amount_bases} and/or --trim_tail2 {amount_bases}


rule fastp_all:
    input:
        expand(str(trim_dir / "{sample}_trimmed.fq.gz"), sample=SAMPLES),
    log: log_dir / "bbduk_all.log"

if config.mode == "pe":
    rule fastp:
        input:
            R1=ancient(str(input_dir / "{sample}.R1.fastq.gz")),
            R2=ancient(str(input_dir / "{sample}.R2.fastq.gz")),
        output:
            R1 = trim_dir / "{sample}.R1.qc.fastq.gz",
            R2 = trim_dir / "{sample}.R2.qc.fastq.gz",
            qhist = log_dir / "{sample}/qhist.txt",
            lhist = log_dir / "{sample}/lhist.txt",
            aqhist = log_dir / "{sample}/aqhist.txt",
        log:
            filter_stats = log_dir / "bbduk.{sample}.log"
        params:
            front_size = 1,
            front_quality = 3,
            tail_size = 1,
            tail_quality = 3,
            right_size = 4,
            right_quality = 15,
            length_required = 36,
        shell:
            'fastp \'
            '--in1 input.R1 \'
            '--in2 input.R2 \'
            '--out1 output.R1 \'
            '--out2 output.R2 \'
            '--json'
            '--html'
            '--cut_front \'
            '--cut_front_window_size params.front_size \'
            '--cut_front_mean_quality params.front_quality \'
            '--cut_tail \'
            '--cut_tail_window_size params.tail_size \'
            '--cut_tail_mean_quality params.tail_quality \'
            '--cut_right \'
            '--cut_right_window_size params.right_size \'
            '--cut_right_mean_quality params.right_quality \'
            '--length_required params.length_required \'

if config.mode == "se":
    rule fastp:
        input:
            R1=ancient(str(input_dir / "{sample}.R1.fastq.gz")),
        output:
            R1 = trim_dir / "{sample}.R1.qc.fastq.gz",
        log:
            filter_stats = log_dir / "bbduk.{sample}.log"
        params:
            front_size = 1,
            front_quality = 3,
            tail_size = 1,
            tail_quality = 3,
            right_size = 4,
            right_quality = 15,
            length_required = 36,
        shell:
            'fastp \'
            '--in1 input.R1 \'
            '--out1 output.R1 \'
            '--json'
            '--html'
            '--cut_front \'
            '--cut_front_window_size params.front_size \'
            '--cut_front_mean_quality params.front_quality \'
            '--cut_tail \'
            '--cut_tail_window_size params.tail_size \'
            '--cut_tail_mean_quality params.tail_quality \'
            '--cut_right \'
            '--cut_right_window_size params.right_size \'
            '--cut_right_mean_quality params.right_quality \'
            '--length_required params.length_required \'


####################################################
### Alignment of reads (paired end & single end) ###
####################################################

rule STAR_make_index:
    input:
        genome = asdf #TODO fix this
    params:
        genome_dir = config.STAR_index,
        annotations = config.annotations,
        overhang = 100,
        threads = config.threads,
        suffix_array_sparsity = 2, # bigger for smaller (RAM), slower indexing
        genomeChrBinNbits = 15, # might need to mess with this for bad genomes
    output:
        directory(params.STAR_index),
    shell:
		'STAR \'
		'--runMode genomeGenerate \'
		'--genomeDir {params.genome_dir} \'
		'--genomeFastaFiles {input.genome} \'
		'--sjdbGTFfile {params.annotations}} \'
		'--sjdbOverhang {params.overhang} \'
		'--runThreadN {params.threads} \'
		'--genomeSAsparseD {params.suffix_array_sparsity} \'
		'--genomeChrBinNbits {params.genomeChrBinNbits}'

rule STAR_all:
    input:
        expand(str(align_dir / "{sample}.aligned.out.bam"), sample=SAMPLES)

if config.mode == "pe":
    rule STAR:
        input:
            R1 = trim_dir / "{sample}.R1.qc.fastq.gz",
            R2 = trim_dir / "{sample}.R2.qc.fastq.gz",
        output:
            bam = align_dir / "{sample}.aligned.out.bam"
        params:
            annotations = config.annotations,
            genome_dir = genome_dir,
            lenfrac_required = config.lenfrac_required
            threads = config.threads,
        resources:
            load=100
        shell:
        	'STAR \'
			'--runThreadN {params.threads} \'
			'--genomeDir {params.genome_dir} \'
			'--readFilesIn {input.R1} {input.R2} \'
			'--quantMode TranscriptomeSAM \'
			'--readFilesCommand zcat \'
			'--outFileNamePrefix {output.bam}' #TODO fix this so that it's pointing to the directory properly

if config.mode == "se":
    rule STAR:
        input:
            R1 = trim_dir / "{sample}.R1.qc.fastq.gz"
        output:
            bam = align_dir / "{sample}.aligned.out.bam"
        params:
            annotations = config.annotations,
            genome_dir = genome_dir,
            lenfrac_required = config.lenfrac_required
            threads = config.threads,
        resources:
            load=100
        shell:
        	'STAR \'
			'--runThreadN {params.threads} \'
			'--genomeDir {params.genome_dir} \'
			'--readFilesIn {input.R1} \'
			'--quantMode TranscriptomeSAM \'
			'--readFilesCommand zcat \'
			'--outFileNamePrefix {output.bam}' #TODO fix this so that it's pointing to the directory properly


#######################
# QUANTIFICATION RSEM #
#######################

rule RSEM_make_index:
    input:
        genome = asdf #TODO fix this
    params:
        index_dir = config.rsem_index_dir,
        annotations = config.annotations,
        genome_id = config.genome_id,
    output:
        directory(params.index_dir),
    shell:
	    'rsem-prepare-reference --gtf {params.annotations} {input.genome} {output} {params.genome_id}'


if config.mode == "pe":
    rule RSEM:
        input:
            bam = align_dir / "{sample}.aligned.out.bam"
        output:
            isoforms = quant_dir / "{sample}.isoforms.results"
            genes = quant_dir / "{sample}.genes.results"
        params:
            threads = config.threads,
        shell:
            'rsem-calculate-expression \'
            '-p {params.threads} \'
            '--paired-end \'
            '--bam {input.bam} \'
            '--no-bam-output \'
            '{input.genome} \'
            '{sample}'

if config.mode == "se":
    rule RSEM:
        input:
            bam = align_dir / "{sample}.aligned.out.bam"
        output:
            isoforms = quant_dir / "{sample}.isoforms.results"
            genes = quant_dir / "{sample}.genes.results"
        params:
            threads = config.threads,
        shell:
            'rsem-calculate-expression \'
            '--paired-end \'
            '--bam {input.bam} \'
            '--no-bam-output \'
            '{input.genome} \'
            '{sample}'

# multiqc

rule multiqc:
    message: "running multiqc"
    input:
        expand(str(fastqc_dir / "{sample}_fastqc.zip"), sample=SAMPLES),
        expand(str(trim_dir / "{sample}_trimmed.fq.gz"), sample=SAMPLES),
        expand(str(align_dir / "{sample}_trimmed_bismark_bt2.bam"), sample=SAMPLES)
    output:
        output_dir / "multiqc_report.html"
    log: log_dir / "multiqc.log"
    run:
        shell("multiqc -fz {input_dir} {trim_dir} {align_dir} {fastqc_dir}")
        shell("mv multiqc_report.html multiqc_data.zip {output_dir}")

