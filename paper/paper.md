---
title: 'bam2fasta: Extract single-cell sequences from bam files to fasta'
authors:
- affiliation: 1
  name: Venkata Naga Pranathi Vemuri
  orcid: 0000-0002-5748-9594
- affiliation: 1
  name: Olga Borisovna Botvinnik
  orcid: 0000-0003-4412-7970

date: "31 December 2019"
output:
  html_document:
    df_print: paged
  pdf_document: default
  word_document: default
bibliography: paper.bib
tags:
- single cell
- bam
- 10x genomics
- fasta
- barcode
affiliations:
- index: 1
  name: Data Sciences Platform, Chan Zuckerberg Biohub, San Francisco, CA
---

# Summary

Single-cell RNA sequencing such as Drop-Seq [McCarroll:2015], 10x Genomics and other microfluidics platforms made leaps over the last decades in the amount of cells that can be sequenced in parallel.
Droplet microfluidics allows Cell Barcodes (CB), along with their unique molecular identifier (UMI) labeled RNA transcripts,  to be sequenced simultaneously with many other cell barcodes of a homogenized tissue.
After alignment, the reads are demultiplexed to cell barcodes, and the whole sequencing run is stored as a binary alignment map file type known as a `.bam` file [Li:2009].
As the demultiplexing occurs after alignment, there is no way to intercept these workflows to extract a single cell's sequences, e.g. as a simple [FASTA])(https://en.wikipedia.org/wiki/FASTA_format) or [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)-formatted file.

As there is no way to identify cells with many reads or few reads *a priori* to alignment, there are many sequences in the bam file with potentially no filter on the minimum number of observations.
There can also be cellular barcodes incorrectly tagged in the `.bam` file due to an error introduced in the chemical reaction, resulting in very few UMI's for that barcode, and these barcodes need to be discarded.
Thus, the size of the resulting `.bam` alignment files are in magnitudes of 10s of GB, and by requiring a minimum number of observed molecular UMIs per cell barcode, could be potentially reduced to 10s of MB, a reduction by three orders of magnitude.
[[Note: would it be possible to add a histogram of the number of reads per UMI? I think that would visually demonstrate the sparsity of data in these files very well.]]

Other existing tools such as `samtools`, `seqtk`, and `bam2bed` currently do not have the ability to remove cellular barcodes with very little data, limiting the downstream potential of the single-cell `.bam` alignment files.
More importantly, the `.bam` may have reads from UMIs that do not contain unique sequence.
Each molecular species is captured and tagged with a UMI, but due to the polymerase chain reaction (PCR) and random shearing of the molecule, different ranges of the sequence could be captured and tagged with the same UMI.
Thus a molecule can appear multiple times with the same UMI, which in other pipelines is less of a problem as they consider only the genomic interval in which the aligned read falls, but in the case of extracting the sequences, then one would want all possible sequences from the same UMI.
Additionally, if there was an error introduced in the molecular processes, there could be several cellular barcodes with the same UMI, causing a barcode collision.
Hence, the erroneous information can be reduced by restricting the data to only the top `n` barcodes with significant number of UMI, which would better represent the genome by removing irregularities and redundancies.
Here we present a tool to filter `.bam` files by CB and UMI and convert them to FASTA files for further data exploration based on sequences per CB.

In this paper we present a technique that converts a binary alignment map `.bam` file to simple sequence FASTA file per cell barcode given threshold such as Unique Molecular Identifier (UMIs) accepted per barcode.`.bam` files can attribute to few limitations as discussed below.
Firstly, loading the UMIs and cell barcodes into memory would require a lot of RAM depending on how the program will allocate memory for differently typed tags in the `.bam` file.
Secondly, recursively iterating through each record in the `.bam` file and merging UMIs within a [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance) is memory intensive, as this method would need a lookup dictionary to be updated as iterations progress, which can case memory leaks and hangups while the huge `.bam` file is still loaded in memory.
Hence we developed the package `bam2fasta` to save the aligned sequence per cell barcode, after the filtering for a mininum UMI, by sharding the `.bam` file into smaller pieces and using Python's `multiprocessing` to clean each piece one by one.

Here is an example of a `.bam` file viewed using the program `samtools` [Li:2009]:

```
$ samtools view $bam | head
A00111:50:H2H5YDMXX:1:1248:13160:2957	16	chr1	138063420	255	1S89M	*	0	0	TTAATAGTTGAAAGTTTATTATGGTTATCAATATTATATCTCAGTAAGAGTAAACAAAACAGTGGGGAAATTCAAGATAAATACACAGTA	F-8FFFFFFFF8FF8FFFFFFFFFFF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFF	NH:i:1	HI:i:1	AS:i:87	nM:i:0	TX:Z:NM_001111316,+5016,89M1S;NM_011210,+4599,89M1S	GX:Z:Ptprc	GN:Z:Ptprc	RE:A:E	CR:Z:AGTTATGTCTCTCGTA	CY:Z:F888-FFFFFF8FFFF	UR:Z:AGGAGGTCTT	UY:Z:F88FFFFFFF	UB:Z:AGGAGGTCTT	BC:Z:CAGTACTG	QT:Z:FFFFFFFF	RG:Z:10X_P7_8:MissingLibrary:1:H2H5YDMXX:1
```

After processing with `bam2fasta` using `--rename-10x-barcodes` flag read name contains the label`lung_epithelial_cell`, CB, and the UMI, this is the resulting FASTA file content:

```
>lung_epithelial_cell|AAATGCCCAAACTGCT-1_GTCATCGCTA_000
TAATAGTTGAAAGTTTATTATGGTTATCAATATTATATCTCAGTAAGAGTAAACAAAACAGTG
```

# Implementation

## Workflow

The `bam2fasta` package solves the issues of memory usage and filtering on UMI content by implementing the following steps.
In the first step, the `.bam` file is sharded into chunks of smaller `.bam` files and is stored in the machine's temporary folder, e.g. `/tmp`.
The chunk size of the `.bam` file is a user-tunable parameter that can be accessed with `--line-count`; by default it is 1500 records.
This process is done serially by iterating through the records in the `.bam` file, using `pysam`, a Python wrapper around `samtools` [Li:2009].

### MapReduce: Map

Fist, we employ a MapReduce [Dean:2008] approach to the temporary `.bam` files to obtain all the reads per cell barcode in a FASTA file.
In the "Map" step, we distribute the computation across multiple processes of the temporary shards of the `.bam` file, such as parsing the barcode, determining the quality of the read, and checking if record is not duplicated.
From these bam shards we create temporary `.fasta` files that contain for each read: the cell barcode, UMI and the aligned sequence.
There might be a cell barcode with a different UMI that would be present in different chunks of these sharded `.bam` files. As a result we would have multiple temporary `.fasta` files for the same barcodes.
We implemented a method to find unique barcodes based on the temporary `.fasta` file names, and for each of the unique barcodes, we assign a temporary barcode `.fasta` files created by different `.bam` shards in a dictionary.

### MapReduce: Reduce

In the "Reduce" step, we combine all sequences for the same barcode. We accomplish this by concatenating strings of temporary FASTA file names for the same barcode, hence its memory consumption is less than it would be if appending to a Python `list` structure.
These temporary `.fasta` files are iteratively split and then combined to one `.fasta` file per barcode by concatenating all the sequences obtained from different `.fasta` files, separated by a user-specified delimiter. 
The default delimiter is `X` as it is not present in common biological alphabets such as DNA or protein alphabets.
For each of the cell barcodes, if their number of UMIs passes the threshold given in flag `--min-umi-per-barcode`, they are considered a "valid" cell barcode.
Using Python's `multiprocessing`, each cellular barcode's FASTAs containing different UMIs are combined into a fasta file with the cellular barcode and the concatenated read sequence, separated by the delimiter `X`, is written to a `.fasta` file.

## Advantages

bam2fasta has several adavantages.
`bam2fasta` can read `.bam` files of any size, and convert to FASTA format quickly.
It is fills the gap to quickly process single-cell RNA-seq `.bam` files, which have unique needs, such as filtering per cell barcode.
This method primarily gives us time and memory performance improvement.
It reduces time from days or just process running out of memory to hours which is concluded from testing on 6-12 GB `.bam` files.
`bam2fasta` takes advantage of sharding which is analogous to tiled rendering in images to save memory and `multiprocessing` and string manipulations to save time.
Depending on the size of `.bam` file and available cores of the compute environment, the time and memory can can be further reduced.


# Installation

The `bam2fasta` package is written in Python, and is available on the [Bioconda](https://bioconda.github.io/) [conda](https://docs.conda.io/en/latest/) channel and the [Python Package Index (PyPI)](https://pypi.org/).
Documentation can be found at https://github.com/czbiohub/bam2fasta/


# Figure

![The bam2fasta workflow as explained in the implementation is illustrated in the flowchart](bam2fasta_workflow.png)


# Glossary

"Sharding," "splitting," "tiling" are synonymously used terms to represent partitioning a complete dataset into smaller subsets.
When the dataset is images, the commonly used term for image rendering world is "tiling."
In computer science, the most common term to explain the phenomenon for any dataset is "sharding."
"Sharding" here is used to enable analyzing a large `.bam` file simultaneously on multiple processes.

"MapReduce" is a phenomenon commonly used in the "Big Data" computing world to map a function/algorithm on a subset of data and reduce/combine the result from each piece of data.

# Acknowledgements

This work was made possible through support from the Chan Zuckerberg Biohub.
Thank you Phoenix Logan, Saba Nafees, and Shayan Hosseinzadeh for helpful discussions.


# References
