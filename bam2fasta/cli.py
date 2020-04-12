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
    if args.method == "shard" and input_format == "bam":
        filenames = tenx_utils.shard_bam_file(
            args.filename,
            args.line_count,
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

        if args.write_barcode_meta_csv:
            tenx_utils.write_to_barcode_meta_csv(
                save_intermediate_files, args.write_barcode_meta_csv)
        fastas = glob.glob(os.path.join(save_fastas, "*_bam2fasta.fasta"))
    else:
        if input_format == "bam":
            basename = os.path.basename(args.filename)
            output_fastq_gzip = os.path.join(
                save_intermediate_files,
                "{}__concatenated.fastq.gz".format(
                    basename.replace(".bam", "")))
            tenx_utils.concatenate_gzip_files(
                [tenx_utils.get_fastq_unaligned(
                    args.filename, n_jobs, save_intermediate_files),
                 tenx_utils.get_fastq_aligned(
                    args.filename, n_jobs, save_intermediate_files)],
                output_fastq_gzip)
        elif input_format == "gz":
            output_fastq_gzip = args.filename

        barcodes_with_significant_umi_records_filename = os.path.join(
            save_intermediate_files,
            "barcodes_with_significant_umi_records.tsv")
        tenx_utils.count_umis_per_cell(
            output_fastq_gzip,
            args.write_barcode_meta_csv,
            args.cell_barcode_pattern,
            args.molecular_barcode_pattern,
            args.min_umi_per_barcode,
            barcodes_with_significant_umi_records_filename)

        unique_barcodes = tenx_utils.read_barcodes_file(
            barcodes_with_significant_umi_records_filename)
        num_unique_barcodes = len(unique_barcodes)
        chunksize = tenx_utils.calculate_chunksize(num_unique_barcodes, n_jobs)
        pool_lists = []
        for i in range(0, num_unique_barcodes, chunksize):
            pool_lists.append(unique_barcodes[i: i + chunksize])
        pool = multiprocessing.Pool(processes=n_jobs)
        logger.info(
            "Pooled %d and chunksize %d mapped for %d lists",
            n_jobs, chunksize, len(pool_lists))

        func = partial(
            tenx_utils.make_per_cell_fastqs,
            output_fastq_gzip,
            args.rename_10x_barcodes,
            save_fastas,
            args.cell_barcode_pattern)

        pool.map(func, pool_lists)

        pool.close()
        pool.join()
        fastas = glob.glob(os.path.join(save_fastas, "*.fastq"))
    logger.info(
        "time taken to write %d fastas is %.5f seconds",
        len(fastas), time.time() - startt)
    return fastas
