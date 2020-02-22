"""
Cli tool to convert 10x bam to fastas
"""

import itertools
import os
import glob
import logging
import time
from collections import OrderedDict
from functools import partial

import screed
from pathos import multiprocessing

from bam2fasta import tenx_utils
from bam2fasta.bam2fasta_args import create_parser, Bam2FastaArgumentParser


logger = logging.getLogger(__name__)

DELIMITER = "X"
CELL_BARCODE = "CELL_BARCODE"
UMI_COUNT = "UMI_COUNT"
READ_COUNT = "READ_COUNT"


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


def info(args):
    "Report bam2fasta version + version of installed dependencies."
    parser = Bam2FastaArgumentParser()
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='report versions of pathos, pysam, and screed')
    args = parser.parse_args(args)

    from bam2fasta import VERSION
    logger.info('bam2fasta version %s', VERSION)
    logger.info('- loaded from path: %s', os.path.dirname(__file__))
    logger.info('')

    if args.verbose:
        import pathos
        logger.info('pathos version %s', pathos.__version__)
        logger.info('- loaded from path: %s', os.path.dirname(pathos.__file__))
        logger.info('')
        import pysam
        logger.info('pysam version %s', pysam.__version__)
        logger.info('- loaded from path: %s', os.path.dirname(pysam.__file__))
        logger.info('')

        logger.info('screed version %s', screed.__version__)
        logger.info('- loaded from path: %s', os.path.dirname(screed.__file__))


def convert(args):
    """Cli tool to convert bam to fasta files"""
    parser = create_parser()
    args = parser.parse_args(args)

    logger.info(args)

    def write_to_barcode_meta_csv():
        """ Merge all the meta text files for each barcode to
        one csv file with CELL_BARCODE, UMI_COUNT,READ_COUNT"""
        barcodes_meta_txts = glob.glob("*_meta.txt")

        with open(args.write_barcode_meta_csv, "w") as fp:
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

    def get_unique_barcodes(all_fastas):
        """ Build a dictionary with each unique barcode as key and
        their fasta files from different shards """
        fasta_files_dict = OrderedDict()
        for fasta in tenx_utils.iter_split(all_fastas, ","):
            barcode = os.path.basename(fasta).replace(".fasta", "")
            value = fasta_files_dict.get(barcode, "")
            fasta_files_dict[barcode] = value + fasta + ","
        # Find unique barcodes
        all_fastas_sorted = list(fasta_files_dict.values())
        all_fastas_sorted.sort()
        del fasta_files_dict
        return all_fastas_sorted

    # Initializing time
    startt = time.time()

    # Setting barcodes file, some 10x files don't have a filtered
    # barcode file
    if args.barcodes_file is not None:
        barcodes = tenx_utils.read_barcodes_file(args.barcodes_file)
    else:
        barcodes = None

    # Shard bam file to smaller bam file
    logger.info('... reading bam file from %s', args.filename)
    n_jobs = args.processes
    filenames = tenx_utils.shard_bam_file(
        args.filename,
        args.line_count,
        args.save_intermediate_files)

    # Create a per-cell fasta generator of sequences
    # If the reads should be filtered by barcodes and umis
    # umis are saved in fasta file as record name and name of
    # the fasta file is the barcode
    func = partial(
        tenx_utils.bam_to_temp_fasta,
        barcodes,
        args.rename_10x_barcodes,
        args.delimiter,
        args.save_intermediate_files)

    length_sharded_bam_files = len(filenames)
    chunksize = calculate_chunksize(length_sharded_bam_files,
                                    n_jobs)
    pool = multiprocessing.Pool(processes=n_jobs)
    logger.info(
        "multiprocessing pool processes {} and chunksize {} calculated".format(
            n_jobs, chunksize))
    # All the fastas are stored in a string instead of a list
    # This saves memory per element of the list by 8 bits
    # If we have unique barcodes in the order of 10^6 before
    # filtering that would result in a huge list if each barcode
    # is saved as a separate element, hence the string
    all_fastas = "," .join(itertools.chain(*(
        pool.imap(
            lambda x: func(x.encode('utf-8')),
            filenames, chunksize=chunksize))))

    # clean up the memmap and sharded intermediary bam files
    [os.unlink(file) for file in filenames if os.path.exists(file)]
    del filenames
    logger.info("Deleted intermediary bam")

    all_fastas_sorted = get_unique_barcodes(all_fastas)
    unique_barcodes = len(all_fastas_sorted)
    logger.info("Found %d unique barcodes", unique_barcodes)
    # Cleaning up to retrieve memory from unused large variables
    del all_fastas
    umi_filter = True if args.min_umi_per_barcode != 0 else False
    if umi_filter:
        tenx_func = tenx_utils.filtered_umi_to_fasta
    else:
        tenx_func = tenx_utils.unfiltered_umi_to_fasta

    chunksize = calculate_chunksize(unique_barcodes, n_jobs)
    pool_lists = []
    for i in range(0, unique_barcodes, chunksize):
        pool_lists.append(all_fastas_sorted[i: i + chunksize])

    logger.info(
        "Pooled %d and chunksize %d mapped for %d lists",
        n_jobs, chunksize, len(pool_lists))

    func = partial(
        tenx_func,
        args.save_fastas,
        args.delimiter,
        args.write_barcode_meta_csv,
        args.min_umi_per_barcode)

    pool.imap(
        lambda pool_list: func(pool_list),
        pool_lists, chunksize=chunksize)

    pool.close()
    pool.join()
    fastas = glob.glob(os.path.join(args.save_fastas, "*_bam2fasta.fasta"))

    if args.write_barcode_meta_csv:
        write_to_barcode_meta_csv()
    logger.info(
        "time taken to convert fastas for 10x folder is %.5f seconds",
        time.time() - startt)
    return fastas
