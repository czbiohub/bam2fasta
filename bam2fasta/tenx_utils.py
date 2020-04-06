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
TENX_TAGS = "CB,UB,XC,XM"

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
    cell_sequences: dict
        dictionary with a cell and corresponding sequence
        ithe cell key is expected to contain umi as well
        separated by the delimiter.
        else {AAAAAAAAAXACTAG: AGCTACACTA} - In this case
        AAAAAAAAA would be cell
        barcode and ACTAG would be umi. The umi will be further
        used by downstream
        processing functions appropriately.
        The barcode is safely returned as the
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
    barcodes_folder = tempfile.mkdtemp(dir=temp_folder)
    for cell, seq in cell_sequences.items():
        barcode, umi = cell.split(delimiter)
        filename = os.path.join(barcodes_folder, barcode + '.fasta')

        # Append to an existing barcode file with a different umi
        with open(filename, "a") as f:
            f.write(">{}\n{}\n".format(umi, seq))
        yield filename


def barcode_umi_seq_to_fasta(
        save_fastas,
        delimiter,
        write_barcode_meta_csv,
        min_umi_per_barcode,
        save_intermediate_files,
        single_barcode_fastas):
    """
    Writes signature records across fasta files for a unique barcode
    Parameters
    ----------
    save_fastas: str
        directory to save the fasta file for the barcode in
    delimiter: str
        separator between two reads, usually 'X'
    write_barcode_meta_csv: bool
        boolean flag, if true
        Metadata per barcode i.e umi count and read count is written
        {barcode}_meta.txt file
    min_umi_per_barcode: int
        Cell barcodes that have less than min_umi_per_barcode
        are ignored
    single_barcode_fastas: str
        comma separated list of fastas belonging to the
        same barcode that were within different bam shards
    save_intermediate_files: str
        Path to save intermediate barcode meta txt files
    """
    # Tracking UMI Counts
    umis = []
    all_split_seqs = []
    # Iterating through fasta files for single barcode from different
    # fastas
    read_count = 0
    logger.info("SINGLE ")
    logger.info(single_barcode_fastas)
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
        unique_meta_file = os.path.join(
            save_intermediate_files, unique_meta_file)
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


def get_fastq_unaligned(input_bam, n_cpus, save_intermediate_files):
    """Get unaligned fastq sequences from bam.

    Parameters
    ----------
    input_bam : str
        Name of the bam file
    n_cpus: int
        number of threads to parallelize the bam file conversion across
    save_intermediate_files: str
        Path to save the output bam and fastq.gz file
    Returns
    -------
    fastq_gz : str
        Path to converted fastq.gz file
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
    """Get aligned fastq sequences from bam.

    Parameters
    ----------
    input_bam : str
        Name of the bam file
    n_cpus: int
        number of threads to parallelize the bam file conversion across
    save_intermediate_files: str
        Path to save the output bam and fastq.gz file
    Returns
    -------
    fastq_gz : str
        Path to converted fastq.gz file
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
    """Return concatenated gzipped file containing
       unaligned and aligned fastq sequences from bam.

    Parameters
    ----------
    input_gz_filenames : list
        list of unaligned and aligned fastq sequences
    Returns
    -------
    output_gz : str
        Path to concatenated fastq.gz file
    """
    with open(output_gz, 'wb') as wfp:
        for fn in input_gz_filenames:
            with open(fn, 'rb') as rfp:
                shutil.copyfileobj(rfp, wfp)


def get_cell_barcode(record, cell_barcode_pattern):
    """Return the cell barcode in the record name.

    Parameters
    ----------
    record : screed record
        screed record containing the cell barcode
    cell_barcode_pattern: regex pattern
        cell barcode pattern to detect in the record name
    Returns
    -------
    barcode : str
        Return cell barcode from the name, if it doesn't exit, returns None
    """
    found_cell_barcode = re.findall(cell_barcode_pattern, record['name'])
    if found_cell_barcode:
        return found_cell_barcode[0][1]


def get_molecular_barcode(record,
                          molecular_barcode_pattern):
    """Return the molecular barcode in the record name.

    Parameters
    ----------
    record : screed record
        screed record containing the molecular barcode
    molecular_barcode_pattern: regex pattern
        molecular barcode pattern to detect in the record name
    Returns
    -------
    barcode : str
        Return molecular barcode from the name,if it doesn't exit, returns None
    """
    found_molecular_barcode = re.findall(molecular_barcode_pattern,
                                         record['name'])
    if found_molecular_barcode:
        return found_molecular_barcode[0][1]


def get_cell_barcode_umis(
        reads,
        cell_barcode_pattern,
        molecular_barcode_pattern):
    """Return a dictionary containing cell barcode string and list of
    corresponding umis as the value

    Parameters
    ----------
    reads : str
        fasta path
    cell_barcode_pattern: regex pattern
        cell barcode pattern to detect in the record name
    molecular_barcode_pattern: regex pattern
        molecular barcode pattern to detect in the record name
    Returns
    -------
    barcode_counter : dict
        dictionary containing cell barcode string and list of
        corresponding umis as the value
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
        barcodes_with_significant_umi_records):
    """Writes to csv the barcodes and number of umis, and to good_barcodes the
    barcodes with greater than or equal to min_umi_per_cell

    Parameters
    ----------
    reads : str
        read records from fasta path
    csv: str
        file to write barcodes and their corresponding number of umis
    cell_barcode_pattern: regex pattern
        cell barcode pattern to detect in the record name
    molecular_barcode_pattern: regex pattern
        molecular barcode pattern to detect in the record name
    min_umi_per_cell: int
        number of minimum umi per cell barcode
    barcodes: str
        csv file containing list of cellular barcode strings
    rename_10x_barcodes : str
        Path to tab-separated file mapping barcodes to their new name
        e.g. with channel or cell annotation label,
        e.g. AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1
    barcodes_with_significant_umi_records: str
        write the valid
        barcodes that have greater than or equal to min_umi_per_cell
    Returns
    -------
    Writes to csv barcodes and their corresponding number of umis
    Writes to barcodes_with_significant_umi_records
    list of the barcodes that have
    greater than or equal to min_umi_per_cell
    """
    barcode_counter = get_cell_barcode_umis(
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
    filtered.to_csv(
        barcodes_with_significant_umi_records, header=False, index=False)


def get_good_cell_barcode_records(
        reads, barcodes_with_significant_umi, cell_barcode_pattern):
    """Returns dict of barcodes as keys and all the records containing
        the barcode with significant umis

    Parameters
    ----------
    reads : str
        fasta path containing all records/reads
    barcodes_with_significant_umi: list
        write the valid
        barcodes that have greater than or equal to min_umi_per_cell
    cell_barcode_pattern: regex pattern
        cell barcode pattern to detect in the record name

    Returns
    -------
    barcodes_with_significant_umi_records: dict
        dict of barcodes as keys and all the records containing
        the barcode with significant umis
    """
    barcodes_with_significant_umi_records = defaultdict(list)

    with screed.open(reads) as f:
        for record in tqdm(f):
            cell_barcode = get_cell_barcode(record, cell_barcode_pattern)
            if cell_barcode in barcodes_with_significant_umi:
                barcodes_with_significant_umi_records[
                    cell_barcode].append(record)
    return barcodes_with_significant_umi_records


def record_to_fastq_string(record):
    """Return the converted fastq string

    Parameters
    ----------
    record : screed record
        record in fasta
    Returns
    -------
    output : str
        convert recotd to fastq string
    """
    return "@{}\n{}\n+\n{}\n".format(
        record['name'], record['sequence'], record['quality'])


def write_fastq(records, filename):
    """Write fastq strings converted records to filename

    Parameters
    ----------
    records : list
        list of screed records
    filename : str
        Path to .fastq string to write to
    Returns
    -------
    Write fastq strings converted records to filename
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


def make_per_cell_fastqs(
        reads,
        barcodes_with_significant_umi_file,
        outdir,
        cell_barcode_pattern,
        num_cpus):
    """Write the filtered cell barcodes in reads
    from barcodes_with_significant_umi_file
    fastq.gzs to outdir

    Parameters
    ----------
    reads : str
        read records from fasta path
    barcodes_with_significant_umi_file: str
        csv file containing barcodes that have
        greater than or equal to min_umi_per_cell
    outdir: str
        write the per cell barcode fastq.gzs to outdir
    cell_barcode_pattern: regex pattern
        cell barcode pattern to detect in the record name
    num_cpus: int
        number of cpus to parallelizing writing per barcode fastas to
    Returns
    -------
    Write the filtered cell barcodes in reads
    from barcodes_with_significant_umi_file
    fastq.gzs to outdir
    """
    barcodes_with_significant_umi = read_barcodes_file(
        barcodes_with_significant_umi_file)

    barcodes_with_significant_umi_records = get_good_cell_barcode_records(
        reads, barcodes_with_significant_umi, cell_barcode_pattern)

    for cell_barcode, records in barcodes_with_significant_umi_records.items():
        filename = os.path.join(outdir, cell_barcode + ".fastq.gz")
        write_fastq(records, filename)
