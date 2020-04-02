"""
10x-sequencing specific utility functions.
"""

import logging
import itertools
import os
from collections import defaultdict
import tempfile
import time

import pysam
import gzip
import shutil
import re
import screed
from tqdm import tqdm
import pandas as pd
import numpy as np

CELL_BARCODES = ['CB', 'XC']
UMIS = ['UB', 'XM']

logger = logging.getLogger(__name__)


def iter_split(string, sep=None):
    """
    Return a generator of strings after
    splitting a string by the given separator

    sep : str
        Separator between strings, default None
    Returns
    -------
    Yields generator of strings after
    splitting a string by the given separator
    """
    sep = sep or ' '
    groups = itertools.groupby(string, lambda s: s != sep)
    return (''.join(g) for k, g in groups if k)


def pass_alignment_qc(alignment, barcodes):
    """
    Check high quality mapping, QC-passing barcode and UMI of alignment.

    alignment :
        aligned bam segment
    barcodes : list
        List of cellular barcode strings
    Returns
    -------
    pass_qc : boolean
        true if a high quality, QC passing barcode with a UMI, false otherwise
    """
    high_quality_mapping = alignment.mapq == 255
    if barcodes is not None:
        good_cell_barcode = any(
            [alignment.has_tag(cb) and alignment.get_tag(cb) in barcodes
             for cb in CELL_BARCODES])
    else:
        good_cell_barcode = any(
            [alignment.has_tag(cb) for cb in CELL_BARCODES])
    good_molecular_barcode = any([alignment.has_tag(umi) for umi in UMIS])
    not_duplicate = not alignment.is_duplicate

    pass_qc = (
        high_quality_mapping and good_cell_barcode and
        good_molecular_barcode and not_duplicate)
    return pass_qc


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


def read_bam_file(bam_path):
    """Read from a QC-pass bam file.

    Parameters
    ----------
    bam_path : str
        Name of a 10x bam file
    Returns
    -------
    bam_file : pysam.AlignmentFile
        Iterator over possorted_genome_bam.bam file
    """
    import pysam

    return pysam.AlignmentFile(bam_path, mode='rb')


def shard_bam_file(bam_file_path, chunked_file_line_count, shards_folder):
    """Shard QC-pass bam file with the given line count
       and save to shards_folder

    Parameters
    ----------
    bam_file_path : str
        Bam file to shard
    chunked_file_line_count: int
        number of lines/alignment reads in each sharded bam file
    shards_folder: str
        absolute path to save the sharded bam files to
    Returns
    -------

    shards : list
        list of sharded bam filenames
    """
    import pysam

    logger.info("Sharding the bam file")
    startt = time.time()
    file_names = []

    with read_bam_file(bam_file_path) as bam_file:
        line_count = 0
        file_count = 0
        header = bam_file.header
        for alignment in tqdm(bam_file):
            if line_count == 0:
                file_name = os.path.join(
                    shards_folder,
                    "temp_bam_shard_{}.bam".format(file_count))
                file_names.append(file_name)
                outf = pysam.AlignmentFile(file_name, "wb", header=header)
            if line_count == chunked_file_line_count:
                file_count = file_count + 1
                line_count = 0
                outf.write(alignment)
                outf.close()
            else:
                outf.write(alignment)
                line_count = line_count + 1
        outf.close()

    logger.info(
        "time taken to shard the bam file into %d shards is %.5f seconds",
        file_count, time.time() - startt)
    return file_names


def bam_to_temp_fasta(
        barcodes, barcode_renamer, delimiter, temp_folder, bam_file):
    """Convert 10x bam to one-record-per-cell fasta.

    Parameters
    ----------
    barcodes : list of str
        QC-passing barcodes
    barcode_renamer : str or None
        Tab-separated filename mapping a barcode to a new name, e.g.
        AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1
    delimiter : str
        Non-DNA or protein alphabet character to be ignored, e.g. if a cell
        has two sequences 'AAAAAAAAA' and 'CCCCCCCC', they would be
        concatenated as 'AAAAAAAAAXCCCCCCCC'.
    temp_folder: str
        folder to save temporary fastas in
    bam : bamnostic.AlignmentFile
    Returns
    -------
    filenames: list
        one temp fasta filename for one cell's high-quality, non-duplicate
        reads

    """
    bam = read_bam_file(bam_file)

    # Filter out high quality alignments and/or alignments with selected
    # barcoddes
    bam_filtered = (x for x in bam if pass_alignment_qc(x, barcodes))
    if barcode_renamer is not None and barcodes is not None:
        renamer = parse_barcode_renamer(barcodes, barcode_renamer)
    else:
        renamer = None
    cell_sequences = defaultdict(str)

    for alignment in bam_filtered:
        # Get barcode of alignment, looks like "AAATGCCCAAACTGCT-1"
        # a bam file might have good cell barcode as any of the tags in
        # CELL_BARCODES
        for cb in CELL_BARCODES:
            if alignment.has_tag(cb):
                barcode = alignment.get_tag(cb)

        renamed = renamer[barcode] if renamer is not None else barcode
        umi = ""
        for umi_tag in UMIS:
            if alignment.has_tag(umi_tag):
                umi = alignment.get_tag(umi_tag)
        renamed = renamed + delimiter + umi

        # Make a long string of all the cell sequences, separated
        # by a non-alphabet letter
        cell_sequences[renamed] += \
            alignment.query_alignment_sequence + delimiter

    filenames = list(set(write_cell_sequences(
        cell_sequences, temp_folder, delimiter)))

    bam.close()

    return filenames


def write_cell_sequences(cell_sequences, temp_folder, delimiter="X"):
    """
    Write each cell's sequences to an individual file
        Parameters
    ----------
    cell_sequences: dictionary with a cell and corresponding sequence
    ithe cell key is expected to contain umi as well
    separated by the delimiter.
    else {AAAAAAAAAXACTAG: AGCTACACTA} - In this case AAAAAAAAA would be cell
    barcode and ACTAG would be umi. The umi will be further used by downstream
    processing functions appropriately. The barcode is safely returned as the
    fasta filename and the umi is saved as record.name/sequence id in the
    fasta file
    delimiter : str, default X
        Used to separate barcode and umi in the cell sequences dict.
    temp_folder: str
        folder to save temporary fastas in

    Returns
    -------
    filenames: generator
        one temp fasta filename for one cell/cell_umi with  sequence
    """
    temp_folder = tempfile.mkdtemp(prefix=temp_folder)
    for cell, seq in cell_sequences.items():
        barcode, umi = cell.split(delimiter)
        filename = os.path.join(temp_folder, barcode + '.fasta')

        # Append to an existing barcode file with a different umi
        with open(filename, "a") as f:
            f.write(">{}\n{}\n".format(umi, seq))
        yield filename


def unfiltered_umi_to_fasta(
        save_fastas,
        delimiter,
        single_barcode_fastas):
    """
    Returns signature records across fasta files for a unique barcode
    Parameters
    ----------
    save_fastas: directory to save the fasta file for the barcode in
    delimiter: separator between two reads, usually 'X'
    single_barcode_fastas: comma separated list of fastas belonging to the
    same barcode that were within different bam shards
    """
    # Getting all fastas for a given barcode
    # from different shards
    count = 0
    # Iterating through fasta files for single barcode from different
    # fastas
    for fasta in iter_split(single_barcode_fastas, ","):

        # Initializing the fasta file to write
        # all the sequences from all bam shards to
        if count == 0:
            unique_fasta_file = os.path.basename(fasta)
            barcode_name = unique_fasta_file.replace(".fasta", "")
            f = open(os.path.join(
                save_fastas, barcode_name + "_bam2fasta.fasta"), "w")

        # Add sequence
        for record in screed.open(fasta):
            sequence = record.sequence
            umi = record.name

            split_seqs = sequence.split(delimiter)
            for index, seq in enumerate(split_seqs):
                if seq == "":
                    continue
                f.write(">{}\n{}\n".format(
                    barcode_name + "_" + umi + "_" + '{:03d}'.format(
                        index), seq))

        # Delete fasta file in tmp folder
        if os.path.exists(fasta):
            os.unlink(fasta)

        count += 1

    # close the fasta file
    f.close()


def filtered_umi_to_fasta(
        save_fastas,
        delimiter,
        write_barcode_meta_csv,
        min_umi_per_barcode,
        single_barcode_fastas):
    """
    Returns signature records across fasta files for a unique barcode
    Parameters
    ----------
    save_fastas: directory to save the fasta file for the barcode in
    delimiter: separator between two reads, usually 'X'
    write_barcode_meta_csv: boolean flag, if true
    Metadata per barcode i.e umi count and read count is written
    {barcode}_meta.txt file
    min_umi_per_barcode: Cell barcodes that have less than min_umi_per_barcode
    are ignored
    single_barcode_fastas: comma separated list of fastas belonging to the
    same barcode that were within different bam shards
    """
    # Tracking UMI Counts
    umis = []
    all_split_seqs = []
    # Iterating through fasta files for single barcode from different
    # fastas
    read_count = 0
    for fasta in iter_split(single_barcode_fastas, ","):
            # calculate unique umi, sequence counts
        for record in screed.open(fasta):
            sequence = record.sequence
            umi = record.name

            # Appending sequence of a umi to the fasta
            split_seqs = sequence.split(delimiter)
            all_split_seqs.append(split_seqs)
            umis.append(umi)
            read_count += len(split_seqs)
        # Delete fasta file in tmp folder
        if os.path.exists(fasta):
            os.unlink(fasta)

    # Write umi count, read count per barcode into a metadata file
    unique_fasta_file = os.path.basename(fasta)
    umi_count = len(umis)
    if write_barcode_meta_csv:
        unique_meta_file = unique_fasta_file.replace(".fasta", "_meta.txt")
        with open(unique_meta_file, "w") as ff:
            ff.write("{} {}".format(umi_count, read_count))

    # If umi count is greater than min_umi_per_barcode write the sequences
    # collected to fasta file for the barcode named as barcode_bam2fasta.fasta
    if umi_count > min_umi_per_barcode:
        barcode_name = unique_fasta_file.replace(".fasta", "")
        with open(
            os.path.join(
                save_fastas,
                barcode_name + "_bam2fasta.fasta"), "w") as f:
            for index, seqs in enumerate(all_split_seqs):
                for seq in seqs:
                    if seqs == "":
                        continue
                    f.write(
                        ">{}\n{}\n".format(
                            barcode_name + "_" +
                            umis[index] + "_" + '{:03d}'.format(index), seq))

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
    basename = os.path.basename(input_bam)
    converted_bam = basename.replace(".bam", "_conveted.bam")
    converted_bam = os.path.join(save_intermediate_files, converted_bam)
    pysam.view(
        input_bam, *["-f4", "-o", converted_bam],
        catch_stdout=False)
    fastq = pysam.fastq(
        converted_bam,
        *["--threads", "{}".format(n_cpus), "-T", TENX_TAGS]).encode()
    fastq_gz = os.path.join(
        save_intermediate_files,
        "{}__unaligned.fastq.gz".format(basename.replace(".bam", "")))
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
    basename = os.path.basename(input_bam)
    converted_bam = basename.replace(".bam", "_conveted.bam")
    converted_bam = os.path.join(save_intermediate_files, converted_bam)
    pysam.view(
        input_bam, *["-ub", "-F", "256", "-q", "255", "-o", converted_bam],
        catch_stdout=False)
    fastq = pysam.fastq(
        converted_bam,
        *["--threads", "{}".format(n_cpus), "-T", TENX_TAGS]).encode()
    fastq_gz = os.path.join(
        save_intermediate_files,
        "{}__aligned.fastq.gz".format(basename.replace(".bam", "")))
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
    if barcodes is not None:
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
