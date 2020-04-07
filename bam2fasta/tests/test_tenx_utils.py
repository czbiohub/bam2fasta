import glob
import itertools
import os
import tempfile

import pandas as pd
import pysam as bs
import pytest
import screed

import bam2fasta.tenx_utils as tenx
import bam2fasta.bam2fasta_args as bam2fasta_args
from bam2fasta.tests import bam2fasta_tst_utils as utils


@pytest.fixture()
def umis():
    return [
        'GCCGTACGGC', 'GGTCGTGGAT', 'ATGTAATAGT', 'ACCGAACGAA', 'CATAACAATA',
        'CCGAGAACCA', 'GGACGGTTTT', 'CATTGCAAGT', 'GGCGCGGCAC', 'GACTAAACTG',
        'TACAACCACG', 'AGGAGGTCTT', 'TGCTAGGAGG', 'GCAAATGGAT', 'CCTAGAACCT',
        'AGCGGCCCAC', 'GCACTCAAGA', 'GACCTTTTAA', 'GTCATCGCTA', 'AATTGACCTG',
        'TTATCACTCG', 'TTAAGACGGG', 'TGGGTATCCT', 'GCGCCAGAGT', 'GATGTTAATT',
        'TGTATCCGGC', 'ACTTCTAGGG', 'CAGTCATTTT', 'CAACCTAGCT', 'GTCAAGTGCT',
        'AGACTATGAA', 'GCACGGAGAC', 'CATAACAATT', 'GGATCGGGAA', 'AATCATGTGG',
        'AGCGGAAATT', 'TGCATCAAGG', 'ACGAGTCCTA', 'GTCGGCAAAT', 'AGAAAATACG',
        'AATGCATGGT', 'GGCCAGCATA', 'AGTAAACAGA', 'GCTGGCCGAT', 'TAGAGAGAGT',
        'AAATGACAAG', 'TACAAATTAA', 'GAACTGGTTG', 'CGGGCAGGGT', 'AAGACTCCTG',
        'AGAATCAATA', 'CCACTTGCAC', 'ATTACAAATG', 'GGGGCGTCTA', 'CCCAAGACGT',
        'AAGACCCAGC', 'ACTCAAACTC', 'ACGGGTTAAG', 'CGTATCCACT', 'AAAGAGTCTT',
        'CGATGTAATG', 'AGGGATCGTG', 'CTAAGTCGCG', 'AGGCACATAT', 'CCACATGCAC',
        'GGCGTAATAC', 'GCTCAGCCCG', 'ACTGTTGACT']


def test_calculate_chunksize():
    tota_jobs_todo = 100
    processes = 1
    obtained = tenx.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = tenx.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained


def test_iter_split():
    expected = ['1', '2', '3']
    input_string = '1,2,3,'
    obtained = list(tenx.iter_split(input_string, ","))
    assert expected == obtained

    obtained = list(tenx.iter_split(input_string, None))
    assert [input_string] == obtained

    input_string = \
        '/path/path2/1.fasta /path/path2/2.fasta /path/path2/3.fasta'
    obtained = list(tenx.iter_split(input_string, None))
    expected = [
        '/path/path2/1.fasta',
        '/path/path2/2.fasta',
        '/path/path2/3.fasta']
    assert expected == obtained


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


def test_barcode_umi_seq_to_fasta_count_0():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    with utils.TempDirectory() as location:
        single_barcode_fastas = tenx.bam_to_temp_fasta(
            barcodes=None,
            barcode_renamer=None,
            delimiter="X",
            bam_file=bam_file,
            temp_folder=location)
        tenx.barcode_umi_seq_to_fasta(
            location,
            "X",
            True,
            0,
            location,
            "," .join(itertools.chain(single_barcode_fastas)))
        fastas = glob.glob(
            os.path.join(location, "*_bam2fasta.fasta"))
        assert len(fastas) == 1


def test_barcode_umi_seq_to_fasta():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    with utils.TempDirectory() as location:
        single_barcode_fastas = tenx.bam_to_temp_fasta(
            barcodes=None,
            barcode_renamer=None,
            delimiter="X",
            bam_file=bam_file,
            temp_folder=location)

        tenx.barcode_umi_seq_to_fasta(
            location,
            "X",
            True,
            10,
            location,
            "," .join(itertools.chain(single_barcode_fastas)))
        fastas = glob.glob(
            os.path.join(location, "*_bam2fasta.fasta"))
        assert len(fastas) == 1
        meta_txts = glob.glob(
            os.path.join(location, "*_meta.txt"))
        assert len(meta_txts) == 1


def test_fastq_unaligned():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    with utils.TempDirectory() as location:
        tenx.get_fastq_unaligned(bam_file, 1, location)
        basename = os.path.basename(bam_file).replace(".bam", "")
        path = os.path.join(
            location, "{}__unaligned.fastq.gz".format(basename))
        assert os.path.exists(path)
        assert os.path.getsize(path) == 58


def test_fastq_aligned():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    with utils.TempDirectory() as location:
        tenx.get_fastq_aligned(bam_file, 1, location)
        basename = os.path.basename(bam_file).replace(".bam", "")
        path = os.path.join(
            location, "{}__aligned.fastq.gz".format(basename))
        assert os.path.exists(path)
        assert os.path.getsize(path) == 50248

        with screed.open(path) as f:
            for record_count, record in enumerate(f):
                assert record != []
        assert record_count == 1708


def test_concatenate_gzip_files():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    with utils.TempDirectory() as location:
        input_gz_filenames = [
            tenx.get_fastq_aligned(bam_file, 1, location),
            tenx.get_fastq_unaligned(bam_file, 1, location)
        ]
        basename = os.path.basename(bam_file).replace(".bam", "")

        path = os.path.join(
            location, "{}__concatenated.fastq.gz".format(basename))
        tenx.concatenate_gzip_files(input_gz_filenames, path)

        assert os.path.exists(path)
        assert os.path.getsize(path) == 50306

        with screed.open(path) as f:
            for record_count, record in enumerate(f):
                assert record != []
        assert record_count == 1708


def test_get_cell_barcode():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    barcodes_file = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = tenx.read_barcodes_file(barcodes_file)
    with utils.TempDirectory() as location:
        tenx.get_fastq_aligned(bam_file, 1, location)
        basename = os.path.basename(bam_file).replace(".bam", "")
        path = os.path.join(
            location, "{}__aligned.fastq.gz".format(basename))
        with screed.open(path) as f:
            for record_count, record in enumerate(f):
                result = tenx.get_cell_barcode(
                    record, bam2fasta_args.CELL_BARCODE_PATTERN)
                if result is not None:
                    assert result in barcodes


def test_get_molecular_barcode(umis):
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    with utils.TempDirectory() as location:
        tenx.get_fastq_aligned(bam_file, 1, location)
        basename = os.path.basename(bam_file).replace(".bam", "")
        path = os.path.join(
            location, "{}__aligned.fastq.gz".format(basename))
        with screed.open(path) as f:
            for record_count, record in enumerate(f):
                result = tenx.get_molecular_barcode(
                    record, bam2fasta_args.MOLECULAR_BARCODE_PATTERN)
                umis.append(result)
                if result is not None:
                    assert result in umis


def test_get_cell_barcode_umis():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    barcodes_file = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = tenx.read_barcodes_file(barcodes_file)
    with utils.TempDirectory() as location:
        tenx.get_fastq_aligned(bam_file, 1, location)
        basename = os.path.basename(bam_file).replace(".bam", "")
        path = os.path.join(
            location, "{}__aligned.fastq.gz".format(basename))
        barcode_cell_umi_dict = tenx.get_cell_barcode_umis(
            path,
            bam2fasta_args.CELL_BARCODE_PATTERN,
            bam2fasta_args.MOLECULAR_BARCODE_PATTERN)
        for barcode in list(barcode_cell_umi_dict.keys()):
            assert barcode in barcodes


def test_count_umis_per_cell():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    barcodes_file = utils.get_test_data('10x-example/barcodes.tsv')
    rename_10x_barcodes = utils.get_test_data(
        '10x-example/barcodes_renamer.tsv')
    with utils.TempDirectory() as location:
        tenx.get_fastq_aligned(bam_file, 1, location)
        basename = os.path.basename(bam_file).replace(".bam", "")
        path = os.path.join(
            location, "{}__aligned.fastq.gz".format(basename))
        meta = os.path.join(location, "barcode_umi_meta.csv")
        good_barcodes = os.path.join(
            location, "barcodes_with_significant_umi_records.csv")
        tenx.count_umis_per_cell(
            path,
            meta,
            bam2fasta_args.CELL_BARCODE_PATTERN,
            bam2fasta_args.MOLECULAR_BARCODE_PATTERN,
            3,
            barcodes_file,
            rename_10x_barcodes,
            good_barcodes)
        df = pd.read_csv(meta)
        df1 = pd.read_csv(good_barcodes)
        print(df.iloc[:, 0].values)
        print(df1.iloc[:, 0].values)
        expected_meta = [15, 2, 2, 5, 4, 6, 2]
        assert expected_meta == pd.read_csv(meta).iloc[:, 0].values.tolist()
        expected_good_barcodes = [
            'lung_epithelial_cell_AAATGCCCAAACTGCT-1',
            'lung_epithelial_cell_AAATGCCGTGAACCTT-1',
            'lung_epithelial_cell_AAATGCCAGATAGTCA-1',
            'human_epithelial_cell_AAACGGGAGGATATAC-1']
        assert expected_good_barcodes == pd.read_csv(
            good_barcodes).iloc[:, 0].values.tolist()


def test_record_to_fastq_string():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    with utils.TempDirectory() as location:
        tenx.get_fastq_aligned(bam_file, 1, location)
        basename = os.path.basename(bam_file).replace(".bam", "")
        path = os.path.join(
            location, "{}__aligned.fastq.gz".format(basename))
        with screed.open(path) as f:
            for record_count, record in enumerate(f):
                result = tenx.record_to_fastq_string(record)
                expected = (
                    "@A00111:50:H2H5YDMXX:2:1334:1886:36072\tUB:Z:AGAAAATACG\n" +
                    "GATTACTTAGTAGCTGTTTACTTAGCAGCACATTTGCAACAGCATCAAAAGCTATGTTACTATAAAATCAGTGCGTGAAGTCTGATTTAC\n")
                assert expected in result
                break


def test_get_good_cell_barcode_records():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    barcodes_file = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = tenx.read_barcodes_file(barcodes_file)
    with utils.TempDirectory() as location:
        tenx.get_fastq_aligned(bam_file, 1, location)
        basename = os.path.basename(bam_file).replace(".bam", "")
        path = os.path.join(
            location, "{}__aligned.fastq.gz".format(basename))
        barcodes_with_significant_umi_records = \
            tenx.get_good_cell_barcode_records(
                path,
                barcodes,
                bam2fasta_args.CELL_BARCODE_PATTERN)
        read_count = 0
        for barcode, records in barcodes_with_significant_umi_records.items():
            read_count += len(records)
            assert barcode in barcodes
        assert read_count == 1610


def test_write_fastq():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    with utils.TempDirectory() as location:
        tenx.get_fastq_aligned(bam_file, 1, location)
        basename = os.path.basename(bam_file).replace(".bam", "")
        path = os.path.join(
            location, "{}__aligned.fastq.gz".format(basename))
        with screed.open(path) as f:
            records = []
            for record_count, record in enumerate(f):
                records.append(record)
        write_path = os.path.join(location, "result.fastq.gz")
        tenx.write_fastq(records, write_path)
        with screed.open(write_path) as f:
            records_written = []
            for record_count, record in enumerate(f):
                records_written.append(record)
            assert records_written == records
        write_path = os.path.join(location, "result.fastq")
        tenx.write_fastq(records, write_path)
        with screed.open(write_path) as f:
            records_written = []
            for record_count, record in enumerate(f):
                records_written.append(record)
            assert records_written == records


def test_make_per_cell_fastqs():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    barcodes_file = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = tenx.read_barcodes_file(barcodes_file)
    with utils.TempDirectory() as location:
        tenx.get_fastq_aligned(bam_file, 1, location)
        basename = os.path.basename(bam_file).replace(".bam", "")
        path = os.path.join(
            location, "{}__aligned.fastq.gz".format(basename))
        outdir = os.path.join(location, "outdir")
        tenx.make_per_cell_fastqs(
            path,
            barcodes_file,
            outdir,
            bam2fasta_args.CELL_BARCODE_PATTERN,
            2)
        fastas = glob.glob(os.path.join(outdir, "*.fastq.gz"))

        for fasta in fastas:
            assert os.path.basename(fasta).replace(".fastq.gz", "") in barcodes
