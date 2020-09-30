bam2fasta
================================
![Tests](https://travis-ci.com/czbiohub/bam2fasta.svg)
[![codecov](https://codecov.io/gh/czbiohub/bam2fasta/branch/master/graph/badge.svg)](https://codecov.io/gh/czbiohub/bam2fasta)

Convert 10x bam file or fastq.gz files to individual FASTA files per cell barcode
convert large bam files to fastq.gz format before the individual fasta files per cell barcode conversion. 
It speeds up this conversion
For small bam files this package can be used directly to convert them to individual FASTA files per cell barcode


To convert bam to fastq.gz format use samtools like below

To get aligned reads from bam file into fastq.gz use 

```
samtools view -ub -F 4 ${bam} \\
    | samtools fastq --threads ${cpus} -T "CB,XC,UB,XM,XB,RG" \\
    | gzip -c - \\
      > ${output_fastq_gz}
```

To get unaligned reads from bam file into fastq.gz use

```
    samtools view -f4 ${bam} \\
      | grep -E '(CB|XC):Z:([ACGT]+)(\\-1)?' \\
      | samtools fastq --threads ${cpus} -T "CB,XC,UB,XM,XB,RG" - \\
      | gzip -c - \\
        > ${output_fastq_gz} \\
```

Using samtools view through python is not recommended for large bam files, as samtools view is streaming the output


Free software: MIT license


## Installation
Latest version can be installed via pip package `bam2fasta`.

Quick install given you have the ssl and zlib packages are already installed.

		pip install bam2fasta
		conda install -c bioconda bam2fasta

Please refer to .travis.yml to see what packages are apt addons on linux and linux addons are required

For osx, before `pip install bam2fasta` install the below homebrew packages

		sudo pip install setuptools
		brew update
		brew install openssl
		brew install zlib

For linux, before `pip install bam2fasta` install the below apt packages

		apt-get install libbz2-dev
		apt-get install libcurl4-openssl-dev
		apt-get install libssl-dev


## Usage

Bam2fasta info command:
	
		bam2fasta info
		bam2fasta info -v

Bam2fasta percell command, it takes BAM and/or barcode files as input. Examples:
	
	bam2fasta percell --filename filename.bam 
	bam2fasta percell --filename 10x-example/possorted_genome_bam.bam \
		--save-fastas fastas --min-umi-per-barcode 10 \
		--write-barcode-meta-csv all_barcodes_meta.csv \
		--barcodes 10x-example/barcodes.tsv \
		--rename-10x-barcodes 10x-example/barcodes_renamer.tsv \
		--shard-size 150 \
		--save-intermediate-files intermediate_files

Bam2fasta percell command, it takes fastq.gz and/or barcode files as input. Examples:
	
	bam2fasta percell --filename 10x-example/possorted_genome_bam.fastq.gz \
		--save-fastas fastas --min-umi-per-barcode 10 \
		--write-barcode-meta-csv all_barcodes_meta.csv \
		--barcodes 10x-example/barcodes.tsv

Bam2fasta count_umis_percell command, it takes fastq.gz file with sequences and barcodes, umis in their read id and counts the umis per cell. Examples:
	
	bam2fasta count_umis_percell --filename filename.fastq.gz 
	bam2fasta count_umis_percell --filename 10x-example/possorted_genome_bam.fastq.gz \
		--write-barcode-meta-csv all_barcodes_meta.csv \
		--min-umi-per-barcode 10 \
		--barcodes-significant-umis-file good_barcodes.csv \
		--cell-barcode-pattern 'CB:Z' \
		--molecular-barcode-pattern 'UB:Z'

Bam2fasta make_fastqs_percell command, it takes it takes fastq.gz file with sequences and barcodes. Examples:
	
	bam2fasta make_fastqs_percell --filename filename.fastq.gz 
	bam2fasta make_fastqs_percell --filename 10x-example/possorted_genome_bam.fastq.gz \
		--save-fastas fastas \
		--min-umi-per-barcode 10 \
		--barcodes-significant-umis-file good_barcodes.csv \
		--barcodes 10x-example/barcodes.tsv \
		--cell-barcode-pattern 'CB:Z'

## Table of Contents
* [Main arguments](#main-arguments)
    * [`--filename`](#--filename)
* [Bam optional parameters](#bam-optional-parameters)
    * [`--barcodes-file`](#--barcodes-file)
    * [`--rename-10x-barcodes`](#--rename-10x-barcodes)
    * [`--save-fastas`](#--save-fastas)
    * [`--save-intermediate-files`](#--save-intermediate-files)
    * [`--write-barcode-meta-csv`](#--write-barcode-meta-csv)
    * [`--min-umi-per-barcode`](#--min-umi-per-barcode)
    * [`--shard-size`](#--shard-size)
    * [`--processes`](#--processes)
    * [`--delimiter`](#--delimiter)
    * [`--cell-barcode-patternt`](#--cell-barcode-pattern)
    * [`--molecular-barcode-pattern`](#--molecular-barcode-pattern)
    * [`--channel-id`](#--channel-id)
    * [`--output-format`](#--output-format)


### `--filename`
For bam/10x files, Use this to specify the location of the bam file or fastq.gz file to get per cell fastas. For example:

```bash
--filename /path/to/data/10x-example/possorted_genome_bam.bam
--filename /path/to/data/10x-example/possorted_genome_bam.fastq.gz
```

## Bam optional parameters


### `--barcodes-file`
For bam/10x files, Use this to specify the location of tsv (tab separated file) containing cell barcodes. For example:

```bash
--barcodes-file /path/to/data/10x-example/barcodes.tsv
```

If left unspecified, barcodes are derived from bam are used.

### `--rename-10x-barcodes`
For bam/10x files, Use this to specify the location of your tsv (tab separated file) containing map of cell barcodes and their corresponding new names(e.g row in the tsv file: AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1). 
For example:

```bash
--rename-10x-barcodes /path/to/data/10x-example/barcodes_renamer.tsv
```
If left unspecified, barcodes in bam as given in barcodes_file are not renamed.


### `--save-fastas`

1. The [save-fastas ](#--save-fastas ) used to save the sequences of each unique barcode in the bam file. By default, they are saved inside working directory to save unique barcodes to files namely {CELL_BARCODE}.fasta. Otherwise absolute path given in save_fastas. 


**Example parameters**

* Default: Save fastas in current working directory:
	* `--save-fastas "fastas"`

### `--save-intermediate-files`

1. The [save-intermediate-files](#--save-intermediate-files ) used to save the intermediate sharded bams and their corresponding fastas. By default, they are saved inside "/tmp/" and are deleted automatically at the end of the program. Otherwise absolute path given in save_intermediate_files. 


**Example parameters**

* Default: Save temporary fastas and bam in `/tmp/` directory:
	* `--save-intermediate-files "fastas"`


### `--write-barcode-meta-csv`
This creates a CSV containing the barcode and number of UMIs per barcode, written in a path given by `write_barcode_meta_csv`. This csv file is empty when the min-umi-per-barcode is zero i.e reads and UMIs per barcode are calculated only when the barcodes are filtered based on [min-umi-per-barcode](#--min-umi-per-barcode)
**Example parameters**

* Default: barcode metadata is not saved 
	* `--write-barcode-meta-csv "barcodes_counts.csv"`


### `--min-umi-per-barcode`
The parameter `--min-umi-per-barcode` ensures that a barcode is only considered a valid barcode read and its sketch is written if number of unique molecular identifiers (UMIs, aka molecular barcodes) are greater than the value specified for a barcode.

**Example parameters**

* Default: min-umi-per-barcode is 0
* Set minimum UMI per cellular barcode as 10:
	* `--min-umi-per-barcode 10`


### `--shard-size`
The parameter `--shard-size` specifies the number of alignments/lines in each bam shard.
**Example parameters**

* Default: shard-size is 1500
	* `--shard-size 400`


### `--processes`
The parameter `--processes` specifies the number of cores/processes to parallelize on.
**Example parameters**

* Default: processes is 2
	* `--processes 400`


### `--delimiter`
The parameter `--delimiter` specifies delimiter between sequences of the same barcode.
**Example parameters**

* Default: delimiter is X
	* `--delimiter X`


### `--cell-barcode-pattern`
The parameter `--cell-barcode-pattern` specifies the regular expressions for molecular barcodes
**Example parameters**

* Default: cell-barcode-pattern is (CB|XC):Z:
	* `--cell-barcode-pattern 'CB:Z'`


### `--molecular-barcode-pattern`
The parameter `--molecular-barcode-pattern` specifies the regular expressions for molecular barcodes.
**Example parameters**

* Default: molecular-barcode-pattern is '(UB|XB):Z:([ACGT]+)'
	* `--molecular-barcode-pattern 'UB:Z'`


### `--channel-id`
The parameter `--channel-id` specifies the prefix for fastqs or fastq.gzs saved by default method.
**Example parameters**

* Default: channel-id is ''
	* `--channel-id 'possorted_aligned'`


### `--output-format`
The parameter `--output-format` specifies the format of output fastq per cell files. it can be either fasta(when input format is bam), fastq or fastq.gz (when input format is fastq.gz). This parameter is only valid for default method
**Example parameters**

* Default: output-format is fastq
	* `--output-format fastq.gz`
