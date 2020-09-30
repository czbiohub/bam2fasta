"""
Cli tool to convert 10x bam to per cell fastas
"""
import itertools
import os
import glob
import logging
import time
from functools import partial

import screed
from pathos import multiprocessing

from bam2fasta import tenx_utils
from bam2fasta.bam2fasta_args import create_parser, Bam2FastaArgumentParser


logger = logging.getLogger(__name__)


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


def count_umis_percell(args):
    if type(args) is list:
        parser = create_parser()
        args = parser.parse_args(args)
        logger.info(args)
    tenx_utils.count_umis_per_cell(
        args.filename,
        args.write_barcode_meta_csv,
        args.cell_barcode_pattern,
        args.molecular_barcode_pattern,
        args.min_umi_per_barcode,
        args.barcodes_significant_umis_file)


def make_fastqs_percell(args):
    if type(args) is list:
        parser = create_parser()
        args = parser.parse_args(args)
        logger.info(args)

    save_path = os.path.abspath(args.save_fastas)
    # Save fasta sequences for aligned and unaligned sequences in
    # separate folders with the following name in save_fastas
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    else:
        logger.info(
            "Path {} already exists, might be overwriting data".format(
                save_path))
    logger.info("Saving fastas at path {}".format(save_path))

    fastqs = tenx_utils.make_per_cell_fastqs(
        args.filename,
        save_path,
        args.channel_id,
        args.output_format,
        args.cell_barcode_pattern,
        args.barcodes_significant_umis_file)

    return fastqs


def percell(args):
    """Cli tool to convert bam to per cell fasta files"""
    parser = create_parser()
    args = parser.parse_args(args)

    logger.info(args)

    save_fastas = os.path.abspath(args.save_fastas)
    if not os.path.exists(save_fastas):
        os.makedirs(save_fastas)
    else:
        logger.info(
            "Path {} already exists, might be overwriting data".format(
                save_fastas))

    # Initializing time
    startt = time.time()

    save_intermediate_files = os.path.abspath(args.save_intermediate_files)
    if not os.path.exists(save_intermediate_files):
        os.makedirs(save_intermediate_files)
    else:
        logger.info(
            "Path {} already exists, might be overwriting data".format(
                save_intermediate_files))

    # Setting barcodes file, some 10x files don't have a filtered
    # barcode file
    if args.barcodes_file is not None:
        barcodes = tenx_utils.read_barcodes_file(args.barcodes_file)
    else:
        barcodes = None

    # Shard bam file to smaller bam file
    logger.info('... reading bam file from %s', args.filename)
    n_jobs = args.processes
    input_format = os.path.basename(args.filename).split(".")[-1]
    if input_format == "bam":
        filenames = tenx_utils.shard_bam_file(
            args.filename,
            args.shard_size,
            save_intermediate_files)

        # Create a per-cell fasta generator of sequences
        # If the reads should be filtered by barcodes and umis
        # umis are saved in fasta file as record name and name of
        # the fasta file is the barcode
        func = partial(
            tenx_utils.bam_to_temp_fasta,
            barcodes,
            args.rename_10x_barcodes,
            args.delimiter,
            save_intermediate_files)

        length_sharded_bam_files = len(filenames)
        chunksize = tenx_utils.calculate_chunksize(
            length_sharded_bam_files,
            n_jobs)
        pool = multiprocessing.Pool(processes=n_jobs)
        logger.info(
            "multiprocessing pool processes {} & chunksize {}".format(
                n_jobs, chunksize))
        # All the fastas are stored in a string instead of a list
        # This saves memory per element of the list by 8 bits
        # If we have unique barcodes in the order of 10^6 before
        # filtering that would result in a huge list if each barcode
        # is saved as a separate element, hence the string
        all_fastas = "," .join(itertools.chain(*(
            pool.imap_unordered(
                lambda x: func(x.encode('utf-8')),
                filenames, chunksize=chunksize))))

        # clean up the memmap and sharded intermediary bam files
        [os.unlink(file) for file in filenames if os.path.exists(file)]
        del filenames
        logger.info("Deleted intermediary bam")

        all_fastas_sorted = tenx_utils.get_fastas_per_unique_barcodes(
            all_fastas)
        unique_barcodes = len(all_fastas_sorted)
        logger.info("Found %d unique barcodes", unique_barcodes)
        # Cleaning up to retrieve memory from unused large variables
        del all_fastas

        # Gather all barcodes per umis to one fasta
        func = partial(
            tenx_utils.barcode_umi_seq_to_fasta,
            save_fastas,
            args.delimiter,
            args.write_barcode_meta_csv,
            args.min_umi_per_barcode,
            save_intermediate_files)

        chunksize = tenx_utils.calculate_chunksize(unique_barcodes, n_jobs)

        logger.info(
            "Pooled %d and chunksize %d mapped for %d lists",
            n_jobs, chunksize, len(all_fastas_sorted))
        list(pool.imap_unordered(
            lambda fasta: func([fasta]),
            all_fastas_sorted, chunksize=chunksize))

        pool.close()
        pool.join()

        # Write barcode meta csv
        if args.write_barcode_meta_csv:
            tenx_utils.write_to_barcode_meta_csv(
                save_intermediate_files, args.write_barcode_meta_csv)

        # Gather all the fastas
        fastas = glob.glob(os.path.join(save_fastas, "*_bam2fasta.fasta"))
    else:
        # if the fastq.gz file is already given
        assert input_format == "gz", \
            ("please convert" +
             "bam files to fastq.gz using samtools")
        # Check if the good barcodes with significant umis is already given
        barcodes_significant_umis_file = args.barcodes_significant_umis_file
        logger.info("barcodes_significant_umis_file {}".format(
            barcodes_significant_umis_file))
        # Find the good_barcodes file from the aligned sequences and use it
        # for unaligned.fastq.gz
        basename_wo_format = os.path.basename(
            args.filename).replace(".fastq.gz", "__")
        args.barcodes_significant_umis_file = os.path.join(
            save_fastas,
            "barcodes_with_significant_umis.tsv")
        count_umis_percell(args)
        args.channel_id = basename_wo_format
        fastas = make_fastqs_percell(args)
        logger.info(
            "time taken to write fastas is %.5f seconds",
            time.time() - startt)

    return fastas


def convert(args):
    # keeping back compatibility
    # bam2fasta 1.0.1 conversion of bam2fsta was called convert
    return percell(args)
