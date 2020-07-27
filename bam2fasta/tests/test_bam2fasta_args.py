import os

from bam2fasta import bam2fasta_args as b2fa_args
from bam2fasta.tests import bam2fasta_tst_utils as utils


def test_bam2fasta_parser():
    parser = b2fa_args.create_parser()
    print(parser)


def test_bam2fasta_valid_args():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example/possorted_genome_bam.bam')
        csv_path = os.path.join(location, "all_barcodes_meta.csv")
        barcodes_path = utils.get_test_data('10x-example/barcodes.tsv')
        renamer_path = utils.get_test_data('10x-example/barcodes_renamer.tsv')
        fastas_dir = os.path.join(location, "fastas")
        save_intermediate_files_dir = os.path.join(location, "temp_fastas")
        if not os.path.exists(fastas_dir):
            os.makedirs(fastas_dir)
        parser = b2fa_args.create_parser()
        args = [
            '--filename', testdata1,
            '--min-umi-per-barcode', '10',
            '--write-barcode-meta-csv', csv_path,
            '--barcodes-file', barcodes_path,
            '--rename-10x-barcodes', renamer_path,
            '--save-fastas', fastas_dir,
            '--save-intermediate-files', save_intermediate_files_dir
        ]
        expected_args_vals = {
            "filename": testdata1,
            "min_umi_per_barcode": 10,
            "write_barcode_meta_csv": csv_path,
            "barcodes_file": barcodes_path,
            "rename_10x_barcodes": renamer_path,
            "save_fastas": fastas_dir,
            "barcodes_significant_umis_file": None,
            "channel_id": "",
            "output_format": "fastq",
            "processes": b2fa_args.DEFAULT_PROCESSES,
            "delimiter": b2fa_args.DEFAULT_DELIMITER,
            "shard_size": b2fa_args.DEFAULT_LINE_COUNT,
            "cell_barcode_pattern": b2fa_args.CELL_BARCODE_PATTERN,
            "molecular_barcode_pattern": b2fa_args.MOLECULAR_BARCODE_PATTERN,
            "save_intermediate_files": save_intermediate_files_dir}
        args = parser.parse_args(args)
        args = vars(args)
        for key, val in args.items():
            assert key in list(expected_args_vals.keys())
            assert val in list(expected_args_vals.values())
