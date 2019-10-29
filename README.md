# bam2fasta
Convert 10x bam file to individual FASTA files of aligned reads per cell

Free software: MIT license


## Installation
Latest version can be installed via pip package `bam2fasta`.

To install this code in developing, clone this github repository and use pip to install

	git clone https://github.com/czbiohub/bam2fasta.git
	cd bam2fasta
	python setup.py develop


## Usage

Both tools have an identical interface and take BAM and/or DataSet files as input. Examples:
	
	bam2fasta filename.bam

* [Main arguments](#main-arguments)
    * [`--filename`](#--filename)
   	* [Bam optional parameters](#bam-optional-parameters)
        * [`--barcodes_file`](#--barcodes_file)
        * [`--rename_10x_barcodes`](#--rename_10x_barcodes)
        * [`--save_fastas`](#--save_fastas)
        * [`--min_umi_per_barcode`](#--min_umi_per_barcode)
        * [`--write_barcode_meta_csv`](#--write_barcode_meta_csv)
        * [`--line_count`](#--line_count)


### `--bam`
For bam/10x files, Use this to specify the location of the bam file. For example:

```bash
--bam /path/to/data/10x-example/possorted_genome_bam
```

## Bam optional parameters


### `--barcodes_file`
For bam/10x files, Use this to specify the location of tsv (tab separated file) containing cell barcodes. For example:

```bash
--barcodes_file /path/to/data/10x-example/barcodes.tsv
```

If left unspecified, barcodes are derived from bam are used.

### `--rename_10x_barcodes`
For bam/10x files, Use this to specify the location of your tsv (tab separated file) containing map of cell barcodes and their corresponding new names(e.g row in the tsv file: AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1). 
For example:

```bash
--rename_10x_barcodes /path/to/data/10x-example/barcodes_renamer.tsv
```
If left unspecified, barcodes in bam as given in barcodes_file are not renamed.


### `--save_fastas`

1. The [save_fastas ](#--save_fastas ) used to save the sequences of each unique barcode in the bam file. It is a path relative to outdir to save unique barcodes to files namely {CELL_BARCODE}.fasta. These fastas are computed once for one permutation of ksize, molecule, and log2_sketch_size, further used to compute the signatures and compare signature matrix for all permutations of ksizes, molecules, and log2_sketch_size. This is done to save the time on saving the computational time and storage in obtaining unique barcodes, sharding the bam file. 


**Example parameters**

* Default: Save fastas in a directory called fastsas inside outdir:
  * `--save_fastas "fastas"`


### `--write_barcode_meta_csv`
This creates a CSV containing the number of reads and number of UMIs per barcode, written in a path relative to `${params.outdir}/barcode_metadata`. This csv file is empty with just header when the min_umi_per_barcode is zero i.e reads and UMIs per barcode are calculated only when the barcodes are filtered based on [min_umi_per_barcode](#--min_umi_per_barcode)
**Example parameters**

* Default: barcode metadata is not saved 
* Save fastas in a file cinside outdir/barcode/metadata:
  * `--write_barcode_meta_csv "barcodes_counts.csv"`


### `--min_umi_per_barcode`
The parameter `--min_umi_per_barcode` ensures that a barcode is only considered a valid barcode read and its sketch is written if number of unique molecular identifiers (UMIs, aka molecular barcodes) are greater than the value specified.

**Example parameters**

* Default: min_umi_per_barcode is 0
* Set minimum UMI per cellular barcode as 10:
  * `--min_umi_per_barcode 10`


### `--line_count`
The parameter `--line_count` specifies the number of alignments/lines in each bam shard.
**Example parameters**

* Default: line_count is 350
* Save fastas in a directory called fastas inside outdir:
  * `--line_count 400`