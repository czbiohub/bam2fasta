import tempfile

import pysam as bs

import bam2fasta.tenx_utils as tenx
from bam2fasta.tests import bam2fasta_tst_utils as utils


def test_read_barcodes_file():
    filename = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = tenx.read_barcodes_file(filename)
    assert len(barcodes) == 10


def test_read_bam_file():
    filename = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    bam_file = tenx.read_bam_file(filename)
    assert isinstance(bam_file, bs.AlignmentFile)
    total_alignments = sum(1 for _ in bam_file)
    assert total_alignments == 1714


def test_shard_bam_file():
    filename = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    bam_file = tenx.read_bam_file(filename)
    assert isinstance(bam_file, bs.AlignmentFile)

    expected_alignments = sum(1 for _ in bam_file)
    with utils.TempDirectory() as location:
        bam_shard_files = tenx.shard_bam_file(
            filename, expected_alignments, location)
        assert len(bam_shard_files) == 1

    num_shards = 2
    with utils.TempDirectory() as location:
        bam_shard_files = tenx.shard_bam_file(
            filename, expected_alignments // num_shards, location)
        assert len(bam_shard_files) == 2

        total_alignments = 0
        for bam_file in bam_shard_files:
            total_alignments += sum(1 for _ in tenx.read_bam_file(bam_file))
        assert total_alignments == expected_alignments

        whole_bam_file = tenx.read_bam_file(filename)
        for bam_file in bam_shard_files:
            for line in tenx.read_bam_file(bam_file):
                assert line == next(whole_bam_file)


def test_pass_alignment_qc():
    barcodes = tenx.read_barcodes_file(
        utils.get_test_data('10x-example/barcodes.tsv'))
    bam = tenx.read_bam_file(
        utils.get_test_data('10x-example/possorted_genome_bam.bam'))

    total_pass = sum(1 for alignment in bam if
                     tenx.pass_alignment_qc(alignment, barcodes))
    assert total_pass == 439


def test_pass_alignment_qc_filtered():
    bam = tenx.read_bam_file(
        utils.get_test_data('10x-example/possorted_genome_bam_filtered.bam'))
    total_alignments = sum(1 for _ in bam)
    bam = tenx.read_bam_file(
        utils.get_test_data('10x-example/possorted_genome_bam_filtered.bam'))
    assert total_alignments == 1500
    total_pass = sum(1 for alignment in bam if
                     tenx.pass_alignment_qc(alignment, None))
    assert total_pass == 192


def test_parse_barcode_renamer():
    filename = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = tenx.read_barcodes_file(filename)
    renamer = tenx.parse_barcode_renamer(barcodes, None)
    for key, value in renamer.items():
        assert key == value
    assert len(renamer) == len(barcodes)

    renamer = tenx.parse_barcode_renamer(
        barcodes, utils.get_test_data('10x-example/barcodes_renamer.tsv'))
    for key, value in renamer.items():
        assert key in value
        assert "epithelial_cell" in value
    assert len(renamer) == len(barcodes)


def test_bam_to_temp_fasta():
    filename = utils.get_test_data('10x-example/barcodes.tsv')
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    barcodes = tenx.read_barcodes_file(filename)

    fastas = tenx.bam_to_temp_fasta(
        barcodes=barcodes,
        barcode_renamer=None,
        delimiter="X",
        bam_file=bam_file,
        temp_folder=tempfile.mkdtemp())
    assert len(list(fastas)) == 8


def test_bam_to_temp_fasta_rename_barcodes():
    filename = utils.get_test_data('10x-example/barcodes.tsv')
    renamer_filename = utils.get_test_data('10x-example/barcodes_renamer.tsv')
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    barcodes = tenx.read_barcodes_file(filename)

    fastas = tenx.bam_to_temp_fasta(
        barcodes=barcodes,
        barcode_renamer=renamer_filename,
        delimiter="X",
        bam_file=bam_file,
        temp_folder=tempfile.mkdtemp())
    assert len(list(fastas)) == 8


def test_filtered_bam_to_umi_fasta():
    bam_file = utils.get_test_data(
        '10x-example/possorted_genome_bam_filtered.bam')
    fastas = tenx.bam_to_temp_fasta(
        barcodes=None,
        barcode_renamer=None,
        delimiter='X',
        bam_file=bam_file,
        temp_folder=tempfile.mkdtemp())
    assert len(list(fastas)) == 32


def test_write_sequences_no_umi():
    cell_sequences = {
        'AAATGCCCAAACTGCT-1X': "atgc",
        'AAATGCCCAAAGTGCT-1X': "gtga"}
    fastas = list(tenx.write_cell_sequences(
        cell_sequences, temp_folder=tempfile.mkdtemp()))
    assert len(fastas) == len(cell_sequences)
    for fasta in fastas:
        assert fasta.endswith(".fasta")


def test_write_sequences_umi():
    cell_sequences = {
        'AAATGCCCAXAACTGCT-1': "atgc",
        'AAATGCCXCAAAGTGCT-1': "gtga",
        'AAATGCCXCAAAGTGCT-2': "gtgc"}
    fastas = list(tenx.write_cell_sequences(
        cell_sequences, temp_folder=tempfile.mkdtemp()))
    assert len(fastas) == len(cell_sequences)
    for fasta in fastas:
        assert fasta.endswith(".fasta")
