# bam2fasta
Convert 10x bam file to individual FASTA files per cell barcode

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

Bam2fasta convert command, it takes BAM and/or barcode files as input. Examples:
	
	bam2fasta percell --filename filename.bam 
	bam2fasta percell --filename 10x-example/possorted_genome_bam.bam \
		--save-fastas fastas --min-umi-per-barcode 10 \
		--write-barcode-meta-csv all_barcodes_meta.csv \
		--barcodes 10x-example/barcodes.tsv \
		--rename-10x-barcodes 10x-example/barcodes_renamer.tsv \
		--line-count 150 \
    --save-intermediate-files intermediate_files

* [Main arguments](#main-arguments)
    * [`--filename`](#--filename)
   	* [Bam optional parameters](#bam-optional-parameters)
        * [`--min-umi-per-barcode`](#--min-umi-per-barcode)
        * [`--write-barcode-meta-csv`](#--write-barcode-meta-csv)
        * [`--processes`](#--processes)
        * [`--delimiter`](#--delimiter)
        * [`--save-fastas`](#--save-fastas)
        * [`--save-intermediate-files`](#--save-intermediate-files)
        * [`--cell-barcode-patternt`](#--cell-barcode-pattern)
        * [`--molecular-barcode-pattern`](#--molecular-barcode-pattern)
        * [`--barcodes-file`](#--barcodes-file)
        * [`--rename-10x-barcodes`](#--rename-10x-barcodes)
        * [`--method`](#--method)


### `--bam`
For bam/10x files, Use this to specify the location of the bam file or tenx.gz file to get per cell fastas/fastqs. For example:

```bash
--bam /path/to/data/10x-example/possorted_genome_bam
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


### `--save_fastas`

1. The [save_fastas ](#--save-fastas ) used to save the sequences of each unique barcode in the bam file. By default, they are saved inside working directory to save unique barcodes to files namely {CELL_BARCODE}.fasta. Otherwise absolute path given in save_fastas. 


**Example parameters**

* Default: Save fastas in current working directory:
  * `--save-fastas "fastas"`

### `--save_intermediate_files`

1. The [save_intermediate_files](#--save-intermediate-files ) used to save the intermediate sharded bams and their corresponding fastas. By default, they are saved inside "/tmp/" and are deleted automatically at the end of the program. Otherwise absolute path given in save_intermediate_files. 


**Example parameters**

* Default: Save temporary fastas and bam in `/tmp/` directory:
  * `--save-intermediate-files "fastas"`


### `--write_barcode_meta_csv`
This creates a CSV containing the number of reads and number of UMIs per barcode, written in a path given by `write_barcode_meta_csv`. This csv file is empty with just header when the min_umi_per_barcode is zero i.e reads and UMIs per barcode are calculated only when the barcodes are filtered based on [min_umi_per_barcode](#--min_umi_per_barcode)
**Example parameters**

* Default: barcode metadata is not saved 
  * `--write-barcode-meta-csv "barcodes_counts.csv"`


### `--min_umi_per_barcode`
The parameter `--min_umi_per_barcode` ensures that a barcode is only considered a valid barcode read and its sketch is written if number of unique molecular identifiers (UMIs, aka molecular barcodes) are greater than the value specified for a barcode.

**Example parameters**

* Default: min_umi_per_barcode is 0
* Set minimum UMI per cellular barcode as 10:
  * `--min-umi-per-barcode 10`


### `--line_count`
The parameter `--line-count` specifies the number of alignments/lines in each bam shard.
**Example parameters**

* Default: line_count is 1500
  * `--line-count 400`


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


### `--cell_barcode_pattern`
The parameter `--cell_barcode_pattern` specifies the regular expressions for molecular barcodes
**Example parameters**

* Default: cell_barcode_pattern is (CB|XC):Z:
  * `--cell-barcode-pattern 'CB:Z'`


### `--molecular_barcode_pattern`
The parameter `--molecular_barcode_pattern` specifies the regular expressions for molecular barcodes.
**Example parameters**

* Default: molecular_barcode_pattern is '(UB|XB):Z:([ACGT]+)'
  * `--molecular-barcode-pattern 'UB:Z'`


### `--method`
The parameter `--method` specifies the method to convert bam to per cell fastas. options: shard bam file and count umis per cell barcode and make per cell fastqs after filtering or do not shard but convert bam to fastq.gz, count umis per cell barcode and filter barcodes and make per cell fastqs.
**Example parameters**

* Default: method is default
  * `--method shard`

