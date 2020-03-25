"""
Cli tool to convert 10x bam to fastas
"""
import os
import logging
import time

import screed

from bam2fasta import tenx_utils
from bam2fasta.bam2fasta_args import create_parser, Bam2FastaArgumentParser


logger = logging.getLogger(__name__)

DELIMITER = "X"
CELL_BARCODE = "CELL_BARCODE"
UMI_COUNT = "UMI_COUNT"
READ_COUNT = "READ_COUNT"


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

    if not os.path.exists(args.save_fastas) and args.save_fastas != "":
        os.makedirs(args.save_fastas)

    # Initializing time
    startt = time.time()

    logger.info('... reading bam file from %s', args.filename)
    n_jobs = args.processes
    save_intermediate_files = args.save_intermediate_files

    output_fastq_gzip = "{}__concatenated.fastq.gz".format(
        args.filename.replace(".bam", ""))
    tenx_utils.concatenate_gzip_files(
        [tenx_utils.get_fastq_unaligned(
            args.filename, n_jobs, save_intermediate_files),
         tenx_utils.get_fastq_aligned(
            args.filename, n_jobs, save_intermediate_files)],
        output_fastq_gzip)

    good_barcodes_filename = os.path.join(
        args.save_intermediate_files, "good_barcodes.tsv")
    tenx_utils.count_umis_per_cell(
        output_fastq_gzip,
        args.write_barcode_meta_csv,
        args.cell_barcode_pattern,
        args.molecular_barcode_pattern,
        args.min_umi_per_barcode,
        args.barcodes_file,
        args.rename_10x_barcodes,
        good_barcodes_filename)

    tenx_utils.make_per_cell_fastas(
        output_fastq_gzip,
        good_barcodes_filename,
        args.save_fastas,
        args.cell_barcode_pattern,
        n_jobs)
