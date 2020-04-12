"""
10x-sequencing specific utility functions.
"""

import logging
import itertools
import os
from collections import defaultdict
import tempfile
import time

import glob
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
CELL_BARCODE = "CELL_BARCODE"
UMI_COUNT = "UMI_COUNT"
READ_COUNT = "READ_COUNT"

logger = logging.getLogger(__name__)


def calculate_chunksize(total_jobs_todo, processes):
    """
    Return a generator of strings after
    splitting a string by the given separator

    sep : str
        Separator between strings, default one space
    Returns
    -------
    Yields generator of strings after
    splitting a string by the given separator
    """
    chunksize, extra = divmod(total_jobs_todo, processes)
    if extra:
        chunksize += 1
    return chunksize


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
    primary = not alignment.is_secondary

    pass_qc = (
        high_quality_mapping and good_cell_barcode and
        good_molecular_barcode and primary)
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
        barcodes = np.unique([line.strip() for line in f]).tolist()
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
        file_count += 1
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
        one temp fasta filename for one cell's high-quality
        reads

    """
    bam = read_bam_file(bam_file)

    # Filter out high quality alignments and/or alignments with selected
    # barcodes
    bam_filtered = (x for x in bam if pass_alignment_qc(x, barcodes))
    if barcode_renamer is not None and barcodes is not None:
        renamer = parse_barcode_renamer(barcodes, barcode_renamer)
    else:
        renamer = None
    cell_sequences = defaultdict(str)

    for count, alignment in enumerate(bam_filtered):
        # Get barcode of alignment, looks like "AAATGCCCAAACTGCT-1"
        # a bam file might have good cell barcode as any of the tags in
        # CELL_BARCODES
        for cb in CELL_BARCODES:
            if alignment.has_tag(cb):
                barcode = alignment.get_tag(cb)
                break

        renamed = renamer[barcode] if renamer is not None else barcode
        umi = ""
        for umi_tag in UMIS:
            if alignment.has_tag(umi_tag):
                umi = alignment.get_tag(umi_tag)
                break
        renamed = renamed + delimiter + umi

        # Make a long string of all the cell sequences, separated
        # by a non-alphabet letter
        cell_sequences[renamed] += \
            alignment.get_forward_sequence() + delimiter
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


def get_fastas_per_unique_barcodes(all_fastas):
    """ Returns the list of fastas per unique barcodes after building
    a dictionary with each unique barcode
    as key and their fasta files from different shards

    Parameters
    ----------
    all_fastas : str
        list of fastas named by barcodes and separated by commas

    Returns
    -------
    Return a list of fastas for all shards per each unique barcode
    """
    fasta_files_dict = defaultdict(str)
    for fasta in iter_split(all_fastas, ","):
        barcode = os.path.basename(fasta).replace(".fasta", "")
        fasta_files_dict[barcode] += fasta + ","
    # Find unique barcodes
    all_fastas_sorted = list(fasta_files_dict.values())
    all_fastas_sorted.sort()
    del fasta_files_dict
    return all_fastas_sorted


def barcode_umi_seq_to_fasta(
        save_fastas,
        delimiter,
        write_barcode_meta_csv,
        min_umi_per_barcode,
        save_files,
        single_barcode_fastas):
    """
    Writes signature records across fasta files for a unique barcode
    Parameters
    ----------
    save_fastas: str
        directory to save the fasta file for the unique barcodes in
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
    save_files: str
        Path to save intermediate barcode meta txt files
    """
    # Tracking UMI Counts
    # Iterating through fasta files for single barcode from different
    # fastas
    for single_barcode_fasta in single_barcode_fastas:
        read_count = 0
        umi_dict = defaultdict(list)
        for fasta in iter_split(single_barcode_fasta, ","):
            # calculate unique umi, sequence counts
            for record in screed.open(fasta):
                sequence = record.sequence
                # Appending sequence of a umi to the fasta
                split_seqs = sequence.split(delimiter)
                if split_seqs[-1] == '':
                    split_seqs = split_seqs[:-1]
                umi_dict[record.name] += split_seqs
                read_count += len(split_seqs)
            # Write umi count, read count per barcode into a metadata file
            unique_fasta_file = os.path.basename(fasta)
            umi_count = len(umi_dict)
            if write_barcode_meta_csv:
                unique_meta_file = unique_fasta_file.replace(
                    ".fasta", "_meta.txt")
                unique_meta_file = os.path.join(
                    save_files, unique_meta_file)
                with open(unique_meta_file, "w") as ff:
                    ff.write("{} {}".format(umi_count, read_count))

            # If umi count is greater than min_umi_per_barcode
            # write the sequences
            # collected to fasta file for barcode named as
            # barcode_bam2fasta.fasta
            # print(fasta, umi_count, read_count, unique_fasta_file)
            if umi_count > min_umi_per_barcode:
                barcode_name = unique_fasta_file.replace(".fasta", "")
                with open(
                    os.path.join(
                        save_fastas,
                        barcode_name + "_bam2fasta.fasta"), "w") as f:
                    for umi, seqs in umi_dict.items():
                        for index, seq in enumerate(seqs):
                            if seq == "":
                                continue
                            f.write(
                                ">{}\n{}\n".format(
                                    barcode_name + "_" +
                                    umi + "_" + '{:03d}'.format(index), seq))


def write_to_barcode_meta_csv(
        barcode_meta_folder, write_barcode_meta_csv):
    """ Merge all the meta text files for each barcode to
    one csv file with CELL_BARCODE, UMI_COUNT,READ_COUNT

    Parameters
    ----------
    barcode_meta_folder : str
        path to folder containing barcode_meta.txt file for all barcodes
        named as barcode.txt and containing umi_count, read_count
    write_barcode_meta_csv : str
        csv file to write the barcode metadata to
    Returns
    -------
    Write csv file to write the barcode metadata to
    """
    barcodes_meta_txts = glob.glob(
        os.path.join(barcode_meta_folder, "*_meta.txt"))
    with open(write_barcode_meta_csv, "w") as fp:
        fp.write("{},{},{}".format(CELL_BARCODE, UMI_COUNT,
                                   READ_COUNT))
        fp.write('\n')
        for barcode_meta_txt in barcodes_meta_txts:
            with open(barcode_meta_txt, 'r') as f:
                umi_count, read_count = f.readline().split()
                umi_count = int(umi_count)
                read_count = int(read_count)

                barcode_name = barcode_meta_txt.replace('_meta.txt', '')
                fp.write("{},{},{}\n".format(barcode_name,
                                             umi_count,
                                             read_count))
            os.unlink(barcode_meta_txt)


def get_fastq_unaligned(input_bam, n_cpus, save_files):
    """Get unaligned fastq sequences from bam.

    Parameters
    ----------
    input_bam : str
        Name of the bam file
    n_cpus: int
        number of threads to parallelize the bam file conversion across
    save_files: str
        Path to save the output bam and fastq.gz file
    Returns
    -------
    fastq_gz : str
        Path to converted fastq.gz file
    """
    basename = os.path.basename(input_bam)
    converted_bam = basename.replace(".bam", "_conveted.bam")
    converted_bam = os.path.join(save_files, converted_bam)
    pysam.view(
        input_bam, *["-f4", "-o", converted_bam],
        catch_stdout=False)
    fastq = pysam.fastq(
        converted_bam,
        *["--threads", "{}".format(n_cpus), "-T", TENX_TAGS]).encode()
    fastq_gz = os.path.join(
        save_files,
        "{}__unaligned.fastq.gz".format(basename.replace(".bam", "")))
    output = gzip.open(fastq_gz, 'wb')
    try:
        output.write(fastq)
    finally:
        output.close()
    return fastq_gz


def get_fastq_aligned(input_bam, n_cpus, save_files):
    """Get aligned fastq sequences from bam.

    Parameters
    ----------
    input_bam : str
        Name of the bam file
    n_cpus: int
        number of threads to parallelize the bam file conversion across
    save_files: str
        Path to save the output bam and fastq.gz file
    Returns
    -------
    fastq_gz : str
        Path to converted fastq.gz file
    """
    basename = os.path.basename(input_bam)
    converted_bam = basename.replace(".bam", "_conveted.bam")
    converted_bam = os.path.join(save_files, converted_bam)
    pysam.view(
        input_bam, *["-ub", "-F", "256", "-q", "255", "-o", converted_bam],
        catch_stdout=False)
    fastq = pysam.fastq(
        converted_bam,
        *["--threads", "{}".format(n_cpus), "-T", TENX_TAGS]).encode()
    fastq_gz = os.path.join(
        save_files,
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
                if molecular_barcode is not None:
                    barcode_counter[cell_barcode].add(molecular_barcode)
    return barcode_counter


def count_umis_per_cell(
        reads,
        csv,
        cell_barcode_pattern,
        molecular_barcode_pattern,
        min_umi_per_cell,
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
        list of the valid
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
                    cell_barcode] += [record]
    return barcodes_with_significant_umi_records


def record_to_fastq_string(record, record_name=None):
    """Return the converted fastq string

    Parameters
    ----------
    record : screed record
        record in fasta
    record_name: str
        specify if you want to rename the record with a new barcode name
    Returns
    -------
    output : str
        convert recotd to fastq string
    """
    if record_name is None:
        result = "@{}\n{}\n+\n{}\n".format(
            record['name'], record['sequence'], record['quality'])
    else:
        result = "@{}\n{}\n+\n{}\n".format(
            record_name, record['sequence'], record['quality'])
    return result


def write_fastq(records, filename, record_name=None):
    """Write fastq strings converted records to filename

    Parameters
    ----------
    records : list
        list of screed records
    filename : str
        Path to .fastq string to write to
    record_name: str
        specify if you want to rename the record with a new barcode name
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
        f.writelines(record_to_fastq_string(r, record_name) for r in records)


def make_per_cell_fastqs(
        reads,
        barcodes_with_significant_umi_file,
        rename_10x_barcodes_file,
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
    rename_10x_barcodes_file: str
        Path to tab-separated file mapping barcodes to their new name
        e.g. with channel or cell annotation label,
        e.g. AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1
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
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        logger.info(
            "Path {} already exists, might be overwriting data".format(outdir))
    barcodes_with_significant_umi = read_barcodes_file(
        barcodes_with_significant_umi_file)

    renamer = parse_barcode_renamer(
        barcodes_with_significant_umi, rename_10x_barcodes_file)

    # Parallelize - subset barcodes and give it to each process
    barcodes_with_significant_umi_records = get_good_cell_barcode_records(
        reads, barcodes_with_significant_umi, cell_barcode_pattern)

    # Parallelize - subset barcodes_with_significant_umi_records and save
    for cell_barcode, records in barcodes_with_significant_umi_records.items():
        renamed = renamer[cell_barcode]
        write_fastq(
            records, os.path.join(outdir, renamed + ".fastq.gz"), renamed)
