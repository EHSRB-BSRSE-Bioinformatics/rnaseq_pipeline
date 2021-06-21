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

# run whole pipeline

rule all:
    input:
        expand(str(trim_dir / "{sample}_trimmed.fq.gz"), sample=SAMPLES),
        expand(str(align_dir / "{sample}_trimmed_bismark_bt2.bam"), sample=SAMPLES),
        output_dir / "multiqc_report.html",
    log: log_dir / "all.log"

# sequence trimming with bbduk

rule bbduk_all:

rule bbduk:


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

