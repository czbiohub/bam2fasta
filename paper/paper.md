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

Next Gen Sequencing like drop-seq, 10x and other microfluidics platforms are making leaps over the last decades in the amount of data it can sequence parallely. With the help of dropseq microfluidics allowing cell barcodes with new unique molecule types to be sequenced simultaneously with many other cell barcodes of a cell type. The whole sequencing run is stored as a `.bam` file. Generally, these workflows don't produce single cell fasta's from the fastq for the whole sequencing run (write why, see kalisto). So, we end up with a lot of data to process in the bam file with no filter on the barcodes and several different UMI's per barcodes. There can also be barcodes incorrectly tagged due to a chemical reaction or not enough quality that enter in the bam file, and they don't have enough UMI's and need to be discarded. Sometimes the size of these `.bam` files tends in the magnitudes of 10s of GB and could be potentially reduced to 10s of MB filtering out the tags to be just the sequence and cell barcode and based on UMI's per barcode. Other existing tools like samtools, seqtk, and bam2bed currently don't have filters to account for that. Therefore we need a tool to filter single cell bam files by cell barcode and convert them to fasta files for further data exploration. In this paper we present a technique that converts a bam file to fasta file per cell barcode given different thresholds like Unique Molecular Identifier (UMIs) accepted per barcode.
Secondly `.bam` files can attribute to few limitations as discussed below. Firstly, loading them in memory all at once would require a lot of RAM depending on how the program will allocate memory for different data typed tags in the `.bam` file. Secondly,recursively going through each record in the `.bam` file and deduce sequence with higher quality and combine sequences with already exisiting barcodes, and different UMIs is memory intensive. This would need a look up dictionary to be updated as it loops through the records in the `.bam` file. This attributes to memory leaks and hangups while the huge `.bam` file is still loaded in memory. Hence the package `bam2fasta` that saves fastas per cell barcode after the filtering using sharding and multiprocessing is implemented.


# Implementation

## Workflow

The package contains solution for the above discussed problems by implementing the following steps.
1. In the first step, sharding the `.bam` file into chunks of smaller `.bam` files and stores them in the machine's temporary folder, e.g. `/tmp`. The chunk size of the `.bam` file is a tunable parameter that can be accessed with `--line_count`; by default it is 1500 records. This process is done serially by iterating through the records in the `.bam` file, using `pysam`, a Python wrapper around samtools [@doi:10.1093/bioinformatics/btp352]. 
### MapReduce: Map

Now we employ a MapReduce [@doi:10.1145/1327452.1327492] approach to the temporary `.bam` files to obtain all the reads per cell barcode in a `.fasta` file.
In the "Map" step, we distribute the computation i.e parsing the barcode, determining the quality of the read, and if record is not duplicated, in parallel across multiple processes on the temporary shards of `.bam` files. These bam shards create temporary `.fasta` files that contain for each read: the cell barcode, UMI and the aligned sequence.
	i.There might be a cell barcode with a different UMI that would be present in different chunks of these sharded `.bam` files. As a result we would have multiple temporary `.fasta` files for the same barcodes. We implemented a method to find unique barcodes based on the temporary `.fasta` file names, and for eacheach of the unique barcodes, assigned a temporary barcode `.fasta` files created by different `.bam` shards in a dictionary.
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


# Acknowledgements

This work was made possible through support from the Chan Zuckerberg Biohub.
Thank you Phoenix Logan (@phoenixAja), Saba Nafees (@snafees), and Shayan Hosseinzadeh (@shayanhoss) for helpful discussions.


# References
