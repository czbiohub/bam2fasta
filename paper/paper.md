---
title: 'bam2fasta: Package to convert bam files to fasta per single cell barcode'
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

Next generation sequencing such as Drop-Seq [@doi:10.1016/j.cell.2015.05.002], 10x Genomics and other microfluidics platforms are making leaps over the last decades in the amount of data it can sequence in parallel. 
Droplet microfluidics allows Cell Barcodes (CB) with their RNA transcripts, labeled by unique molecule identifiers (UMIs) to be sequenced simultaneously with many other cell barcodes of a homogenized tissue. 
After alignment, the cell barcodes are demultiplexed and whole sequencing run is stored as a binary alignment map file type known as a `.bam` file [@doi:10.1093/bioinformatics/btp352]. 
As the demultiplexing occurs after alignment in every single-cell RNA-seq pipeline, these workflows don't produce single cell fasta's from the fastq for the whole sequencing run. 

Thus, there are many RNA sequences in the bam file with potentially no filter on the minimum number of observations.
There can also be cellular barcodes incorrectly tagged due to a chemical reaction that are present in the bam file, and as they do not have enough UMI's and need to be discarded.
This results in the size of the `.bam` files in the magnitudes of 10s of GB and could be potentially reduced to 10s of MB by requiring a minimum number of observed molecular UMIs per cell barcode. 

Other existing tools like `samtools`, `seqtk`, and `bam2bed` currently don't have filters to account for the fact that most of the tags depending on further study are not essential. 
More importantly the tools have several CBs and UMIs that are not unique. 
Each UMI or barcode can be captured and tagged multiple times and can be entered as a read, some genes which do not have many UMIs can be tagged incorrectly and saved as a read with barcode. 
There could be several barcodes with single UMIs if there is an error. 
Hence this information can be reduced to the top `n` barcodes with significant number of UMI. 
This could better represent the genome by removing irregularities and redundancies. 
Hence we present a tool to filter bam files by CB and UMI and convert them to fasta files for further data exploration based on sequences per CB. 

In this paper we present a technique that converts a binary alignment map `.bam` file to simple sequence `.fasta` file per cell barcode given threshold like Unique Molecular Identifier (UMIs) accepted per barcode.`.bam` files can attribute to few limitations as discussed below. 
Firstly, loading them in memory all at once would require a lot of RAM depending on how the program will allocate memory for different data typed tags in the `.bam` file. 
Secondly,recursively going through each record in the `.bam` file and deduce sequence with higher quality and combine sequences with already exisiting barcodes, and different UMIs is memory intensive. This would need a look up dictionary to be updated as it loops through the records in the `.bam` file. This attributes to memory leaks and hangups while the huge `.bam` file is still loaded in memory. 
Hence the package `bam2fasta` that saves fastas per cell barcode after the filtering using sharding and multiprocessing is implemented.

Here is an example of the `.bam` file and using samtools
`samtools view $bam | head`
`A00111:50:H2H5YDMXX:1:1248:13160:2957	16	chr1	138063420	255	1S89M	*	0	0	TTAATAGTTGAAAGTTTATTATGGTTATCAATATTATATCTCAGTAAGAGTAAACAAAACAGTGGGGAAATTCAAGATAAATACACAGTA	F-8FFFFFFFF8FF8FFFFFFFFFFF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFF	NH:i:1	HI:i:1	AS:i:87	nM:i:0	TX:Z:NM_001111316,+5016,89M1S;NM_011210,+4599,89M1S	GX:Z:Ptprc	GN:Z:Ptprc	RE:A:E	CR:Z:AGTTATGTCTCTCGTA	CY:Z:F888-FFFFFF8FFFF	UR:Z:AGGAGGTCTT	UY:Z:F88FFFFFFF	UB:Z:AGGAGGTCTT	BC:Z:CAGTACTG	QT:Z:FFFFFFFF	RG:Z:10X_P7_8:MissingLibrary:1:H2H5YDMXX:1`

The resulting `.fasta` would look like
`>lung_epithelial_cell|AAATGCCCAAACTGCT-1_GTCATCGCTA_000
TAATAGTTGAAAGTTTATTATGGTTATCAATATTATATCTCAGTAAGAGTAAACAAAACAGTG`

# Implementation

## Workflow

The package contains solution for the above discussed problems by implementing the following steps.
In the first step, sharding the `.bam` file into chunks of smaller `.bam` files and stores them in the machine's temporary folder, e.g. `/tmp`. 
The chunk size of the `.bam` file is a tunable parameter that can be accessed with `--line-count`; by default it is 1500 records. 
This process is done serially by iterating through the records in the `.bam` file, using `pysam`, a Python wrapper around samtools [@doi:10.1093/bioinformatics/btp352]. 

### MapReduce: Map

Now we employ a MapReduce [@doi:10.1145/1327452.1327492] approach to the temporary `.bam` files to obtain all the reads per cell barcode in a `.fasta` file.
In the "Map" step, we distribute the computation i.e parsing the barcode, determining the quality of the read, and if record is not duplicated, in parallel across multiple processes on the temporary shards of `.bam` files. These bam shards create temporary `.fasta` files that contain for each read: the cell barcode, UMI and the aligned sequence.
There might be a cell barcode with a different UMI that would be present in different chunks of these sharded `.bam` files. As a result we would have multiple temporary `.fasta` files for the same barcodes. 
We implemented a method to find unique barcodes based on the temporary `.fasta` file names, and for each of the unique barcodes, assigned a temporary barcode `.fasta` files created by different `.bam` shards in a dictionary.

### MapReduce: Reduce

In the "Reduce" step, we combine all sequences for the same barcode.

We accomplish this by concatenating strings of temporary `.fasta` file names for the same barcode, hence its memory consumption is less than it would be if appending to a Python `list` structure. 
These temporary `.fasta` files are iteratively split and then combined to one `.fasta` file per barcode by concatenating all the sequences obtained from different `.fasta` files, separated by a user-specified delimiter. 
For each of the cell barcodes, only valid cell barcodes are obtained based on the number of UMI's per cell barcode given in flag `--min-umi-per-barcode`. 
Each barcodes with a subset of `.fasta`s with different UMI's are combined simultaneously and a fasta file with barcode and the concatenated sequence separated by `X` is written to a `.fasta` file using multiprocessing.

## Advantages

bam2fasta has several adavantages.
`bam2fasta` can read `.bam` files of any size, and convert to FASTA format quickly. 
It is fills the gap to quickly process single-cell RNA-seq `.bam` files, which have unique needs, such as filtering per cell barcode.
This method primarily gives us time and memory performance improvement. 
It reduces time from days or just process running out of memory to hours which is concluded from testing on 6-12 GB `.bam` files. 
`bam2fasta` takes advantage of sharding which is analogous to tiled rendering in images to save memory and multiprocessing and string manipulations to save time. 
Depending on the size of `.bam` file and resources of the cluster/computer it can be further reduced.


# Installation

The `bam2fasta` package is written in Python, and is available on the [Bioconda](https://bioconda.github.io/) [conda](https://docs.conda.io/en/latest/) channel and the [Python Package Index (PyPI)](https://pypi.org/).
Documentation can be found at https://github.com/czbiohub/bam2fasta/


# Figure

![The bam2fasta workflow as explained in the implementation is illustrated in the flowchart](bam2fasta_workflow.png)


# Glossary

Sharding, splitting, tiling synonymously used terms represent looking at a subset of a complete dataset. When the dataset is images it's usually referred in image rendering world as tiling. In computer science, the most common term to explain the phenomenon for any dataset is sharding. 
Sharding here is used to enable analyzing a large `.bam` file simultaneously on multiple processes.

MapReduce is a phenomenon commonly used in the "Big Data" computing world to map a function/algorithm on a subset of data and reduce/combine the result from each piece of data.

# Acknowledgements

This work was made possible through support from the Chan Zuckerberg Biohub.
Thank you Phoenix Logan (@phoenixAja), Saba Nafees (@snafees), and Shayan Hosseinzadeh (@shayanhoss) for helpful discussions.


# References
