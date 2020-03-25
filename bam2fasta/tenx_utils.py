"""
10x-sequencing specific utility functions.
"""

import logging
from collections import defaultdict
import sys
import itertools
import multiprocessing
import time
import re
import os
from io import BufferedReader

import gzip
import shutil
import screed
import pysam
import pandas as pd
from tqdm import tqdm
import numpy as np

TENX_TAGS = "CB,UB,XC,XM"
logger = logging.getLogger(__name__)


def get_fastq_unaligned(input_bam, n_cpus, save_intermediate_files):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    converted_bam = input_bam.replace(".bam", "_conveted.bam")
    pysam.view(
        input_bam, *["-f4", converted_bam],
        catch_stdout=False)
    fastq = pysam.fastq(
        converted_bam,
        *["--threads", "{}".format(n_cpus), "-T", TENX_TAGS]).encode()
    fastq_gz = os.path.join(
        save_intermediate_files,
        "{}__unaligned.fastq.gz".format(input_bam.replace(".bam", "")))
    output = gzip.open(fastq_gz, 'wb')
    try:
        output.write(fastq)
    finally:
        output.close()
    return fastq_gz


def get_fastq_aligned(input_bam, n_cpus, save_intermediate_files):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    converted_bam = input_bam.replace(".bam", "_conveted.bam")
    pysam.view(
        input_bam, *["-ub", "-F", "256", "-q", "255", "-o", converted_bam],
        catch_stdout=False)
    fastq = pysam.fastq(
        converted_bam,
        *["--threads", "{}".format(n_cpus), "-T", TENX_TAGS]).encode()
    fastq_gz = os.path.join(
        save_intermediate_files,
        "{}__aligned.fastq.gz".format(input_bam.replace(".bam", "")))

    output = gzip.open(fastq_gz, 'wb')
    try:
        output.write(fastq)
    finally:
        output.close()
    return fastq_gz


def concatenate_gzip_files(input_gz_filenames, output_gz):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    with open(output_gz, 'wb') as wfp:
        for fn in input_gz_filenames:
            with open(fn, 'rb') as rfp:
                shutil.copyfileobj(rfp, wfp)


def get_cell_barcode(record, cell_barcode_pattern):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    found_cell_barcode = re.findall(cell_barcode_pattern, record['name'])
    if found_cell_barcode:
        return found_cell_barcode[0][1]


def get_molecular_barcode(record,
                          molecular_barcode_pattern):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    found_molecular_barcode = re.findall(molecular_barcode_pattern,
                                         record['name'])
    if found_molecular_barcode:
        return found_molecular_barcode[0][1]


def get_cell_barcode_umi_counts(reads,
                                cell_barcode_pattern,
                                molecular_barcode_pattern):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    barcode_counter = defaultdict(set)

    with screed.open(reads) as f:
        for record in tqdm(f):
            cell_barcode = get_cell_barcode(record, cell_barcode_pattern)
            if cell_barcode is not None:
                molecular_barcode = get_molecular_barcode(
                    record,
                    molecular_barcode_pattern)
                barcode_counter[cell_barcode].add(molecular_barcode)
    return barcode_counter


def parse_barcode_renamer(barcodes, barcode_renamer):
    """
    Return a dictionary with cell barcode and the renamed barcode.

    barcodes : list
        List of cellular barcode strings
    barcode_renamer : str
        Path to tab-separated file mapping barcodes to their new name
        e.g. with channel or cell annotation label,
        e.g. AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1
    Returns
    -------
    barcode_renamer : dict
        A (str, str) mapping of the original barcode to its new name
    """
    if barcode_renamer is not None:
        renamer = {}

        with open(barcode_renamer) as f:
            for line in f.readlines():
                barcode, renamed = line.split()
                assert barcode in barcodes
                renamer[barcode] = renamed.replace("|", "_")
    else:
        renamer = dict(zip(barcodes, barcodes))
    return renamer


def read_barcodes_file(barcode_path):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    with open(barcode_path) as f:
        barcodes = np.unique([line.strip() for line in f])
    return barcodes


def count_umis_per_cell(
        reads,
        csv,
        cell_barcode_pattern,
        molecular_barcode_pattern,
        min_umi_per_cell,
        barcodes,
        rename_10x_barcodes,
        good_barcodes):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    barcode_counter = get_cell_barcode_umi_counts(
        reads,
        cell_barcode_pattern,
        molecular_barcode_pattern)
    if barcodes is None:
        renamer = parse_barcode_renamer(
            read_barcodes_file(barcodes), rename_10x_barcodes)
        umi_per_barcode = {
            renamer[k]: len(v) for k, v in barcode_counter.items()}
    else:
        umi_per_barcode = {
            k: len(v) for k, v in barcode_counter.items()}
    series = pd.Series(umi_per_barcode)
    series.to_csv(csv, header=False, index=False)

    filtered = pd.Series(series[series >= min_umi_per_cell].index)
    filtered.to_csv(good_barcodes, header=False, index=False)


def get_good_cell_barcode_records(reads, good_barcodes, cell_barcode_pattern):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    good_cell_barcode_records = defaultdict(list)

    with screed.open(reads) as f:
        for record in tqdm(f):
            cell_barcode = get_cell_barcode(record, cell_barcode_pattern)
            if cell_barcode in good_barcodes:
                good_cell_barcode_records[cell_barcode].append(record)
    return good_cell_barcode_records


def record_to_fastq_string(record):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    return "@{}\n{}\n+\n{}\n".format(
        record['name'], record['sequence'], record['quality'])


def write_records(records, filename):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    if filename.endswith('gz'):
        import gzip
        opener = gzip.open
        mode = 'wt'
    else:
        opener = open
        mode = 'w'

    with opener(filename, mode) as f:
        f.writelines([record_to_fastq_string(r) for r in records])


def make_per_cell_fastas(
        reads,
        good_barcodes_filename,
        outdir,
        cell_barcode_pattern,
        num_cpus):
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    good_barcodes = read_barcodes_file(good_barcodes_filename)

    good_cell_barcode_records = get_good_cell_barcode_records(
        reads, good_barcodes, cell_barcode_pattern)

    for cell_barcode, records in good_cell_barcode_records.items():
        filename = os.path.join(outdir, cell_barcode + ".fastq.gz")
        write_records(records, filename)
