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

# sequence trimming with bbduk

rule bbduk_all:
    input:
        expand(str(trim_dir / "{sample}_trimmed.fq.gz"), sample=SAMPLES),
    log: log_dir / "bbduk_all.log"

if config.mode == "pe":
    rule bbduk:
        input:
            R1=ancient(str(input_dir / "{sample}.R1.fastq.gz")),
            R2=ancient(str(input_dir / "{sample}.R2.fastq.gz")),
        output:
            R1=trim_dir / "{sample}.R1.qc.fastq.gz" %(indir),
            R2=trim_dir / "{sample}.R2.qc.fastq.gz" %(indir),
            qhist= log_dir / "{sample}/qhist.txt" %(indir),
            lhist= log_dir / "{sample}/lhist.txt" %(indir),
            aqhist= log_dir / "{sample}/aqhist.txt" %(indir),
        log:
            filter_stats = log_dir / "bbduk.{sample}.log"
        params:
            qtrim=config.qtrim,
            quality=config.quality,
            min_len=config.min_len,
            adaptors=config.adaptors,
        shell:
            '''
            bbduk.sh in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} maq={params.quality} ref={params.adaptors} qtrim={params.qtrim} trimq={params.quality} minlength={params.min_len} tpe tbo qhist={output.qhist} lhist={output.lhist} aqhist={output.aqhist} overwrite=t 2> {log.filter_stats}
            '''
            # If you put this into a run: with a shell() it will complain about exitcode != 0 when you use 2> {log}
            # See https://snakemake.readthedocs.io/en/stable/project_info/faq.html#my-shell-command-fails-with-exit-code-0-from-within-a-pipe-what-s-wrong

if config.mode == "se":
    rule bbduk:
        input:
            R1=ancient(str(input_dir / "{sample}.R1.fastq.gz")),
        output:
            R1=trim_dir / "{sample}.R1.qc.fastq.gz" %(indir),
            qhist= log_dir / "{sample}/qhist.txt" %(indir),
            lhist= log_dir / "{sample}/lhist.txt" %(indir),
            aqhist= log_dir / "{sample}/aqhist.txt" %(indir),
        log:
            filter_stats = log_dir / "bbduk.{sample}.log"
        params:
            qtrim=config.qtrim,
            quality=config.quality,
            min_len=config.min_len,
            adaptors=config.adaptors,
        shell:
            '''
            bbduk.sh in={input.R1} out={output.R1} maq={params.quality} ref={params.adaptors} qtrim={params.qtrim} trimq={params.quality} minlength={params.min_len} tpe tbo qhist={output.qhist} lhist={output.lhist} aqhist={output.aqhist} overwrite=t 2> {log.filter_stats}
            '''
            # If you put this into a run: with a shell() it will complain about exitcode != 0 when you use 2> {log}
            # See https://snakemake.readthedocs.io/en/stable/project_info/faq.html#my-shell-command-fails-with-exit-code-0-from-within-a-pipe-what-s-wrong


# alignment with STAR

rule STAR_all:

rule STAR:

# fastqc

rule fastqc_all:
    input:
        expand(str(fastqc_dir / "{sample}_fastqc.zip"), sample=SAMPLES)

rule fastqc:
    input:
        input_dir / "{sample}.fastq.gz"
    output: 
        fastqc_dir / "{sample}_fastqc.zip"
    run:
        shell("fastqc -o {fastqc_dir} {input}")


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

