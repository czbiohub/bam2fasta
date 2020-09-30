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


@pytest.fixture()
def expected_good_barcodes():
    return [
        'AAAGATGCAGATCTGT-1',
        'AAATGCCCAAACTGCT-1',
        'AAATGCCGTGAACCTT-1',
        'AAATGCCAGATAGTCA-1',
        'AAACGGGAGGATATAC-1']


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
    assert total_pass == 1610


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


def test_get_fastas_per_unique_barcodes():
    filename = utils.get_test_data('10x-example/barcodes.tsv')
    renamer_filename = utils.get_test_data('10x-example/barcodes_renamer.tsv')
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    barcodes = tenx.read_barcodes_file(filename)

    with utils.TempDirectory() as location:

        all_fastas = tenx.bam_to_temp_fasta(
            barcodes=barcodes,
            barcode_renamer=renamer_filename,
            delimiter="X",
            bam_file=bam_file,
            temp_folder=location)
        all_fastas = ",".join(itertools.chain(all_fastas))
        fastas_sorted = tenx.get_fastas_per_unique_barcodes(all_fastas)
        assert len(fastas_sorted) == 8


def test_barcode_umi_seq_to_fasta_count_0():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    with utils.TempDirectory() as location:
        single_barcode_fastas = tenx.bam_to_temp_fasta(
            barcodes=None,
            barcode_renamer=None,
            delimiter="X",
            bam_file=bam_file,
            temp_folder=location)
        single_barcode_fastas = "," .join(
            itertools.chain(single_barcode_fastas))
        all_fastas_sorted = tenx.get_fastas_per_unique_barcodes(
            single_barcode_fastas)
        tenx.barcode_umi_seq_to_fasta(
            location,
            "X",
            True,
            0,
            location,
            all_fastas_sorted)
        fastas = glob.glob(
            os.path.join(location, "*_bam2fasta.fasta"))
        assert len(fastas) == 8


def test_write_to_barcode_meta_csv():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    with utils.TempDirectory() as location:
        single_barcode_fastas = tenx.bam_to_temp_fasta(
            barcodes=None,
            barcode_renamer=None,
            delimiter="X",
            bam_file=bam_file,
            temp_folder=location)
        single_barcode_fastas = "," .join(
            itertools.chain(single_barcode_fastas))
        all_fastas_sorted = tenx.get_fastas_per_unique_barcodes(
            single_barcode_fastas)
        tenx.barcode_umi_seq_to_fasta(
            location,
            "X",
            True,
            0,
            location,
            all_fastas_sorted)
        csv = os.path.join(location, "meta.csv")
        tenx.write_to_barcode_meta_csv(location, csv)
        umi_counts = [6, 2, 6, 4, 15, 5, 2, 2]
        read_counts = [312, 36, 251, 194, 594, 153, 2, 68]
        for index, row in pd.read_csv(csv).iterrows():
            assert umi_counts[index] == row[tenx.UMI_COUNT]
            assert read_counts[index] == row[tenx.READ_COUNT]


def test_barcode_umi_seq_to_fasta():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    with utils.TempDirectory() as location:
        single_barcode_fastas = tenx.bam_to_temp_fasta(
            barcodes=None,
            barcode_renamer=None,
            delimiter="X",
            bam_file=bam_file,
            temp_folder=location)
        single_barcode_fastas = "," .join(
            itertools.chain(single_barcode_fastas))
        all_fastas_sorted = tenx.get_fastas_per_unique_barcodes(
            single_barcode_fastas)
        tenx.barcode_umi_seq_to_fasta(
            location,
            "X",
            True,
            10,
            location,
            all_fastas_sorted)
        fastas = glob.glob(
            os.path.join(location, "*_bam2fasta.fasta"))
        assert len(fastas) == 1
        meta_txts = glob.glob(
            os.path.join(location, "*_meta.txt"))
        assert len(meta_txts) == 8


def test_get_cell_barcode():
    fastq_file = utils.get_test_data(
        '10x-example/possorted_genome_bam.fastq.gz')
    barcodes_file = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = tenx.read_barcodes_file(barcodes_file)
    with screed.open(fastq_file) as f:
        for record_count, record in enumerate(f):
            result = tenx.get_cell_barcode(
                record, bam2fasta_args.CELL_BARCODE_PATTERN)
            if result is not None:
                assert result in barcodes


def test_get_molecular_barcode(umis):
    fastq_file = utils.get_test_data(
        '10x-example/possorted_genome_bam.fastq.gz')
    with screed.open(fastq_file) as f:
        for record_count, record in enumerate(f):
            result = tenx.get_molecular_barcode(
                record, bam2fasta_args.MOLECULAR_BARCODE_PATTERN)
            umis.append(result)
            if result is not None:
                assert result in umis


def test_get_cell_barcode_umis():
    path = utils.get_test_data('10x-example/possorted_genome_bam.fastq.gz')
    barcodes_file = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = tenx.read_barcodes_file(barcodes_file)
    barcode_cell_umi_dict = tenx.get_cell_barcode_umis(
        path,
        bam2fasta_args.CELL_BARCODE_PATTERN,
        bam2fasta_args.MOLECULAR_BARCODE_PATTERN)
    for barcode in list(barcode_cell_umi_dict.keys()):
        assert barcode in barcodes


def test_count_umis_per_cell(expected_good_barcodes):
    path = utils.get_test_data('10x-example/possorted_genome_bam.fastq.gz')
    with utils.TempDirectory() as location:
        meta = os.path.join(location, "barcode_umi_meta.csv")
        good_barcodes = os.path.join(
            location, "barcodes_with_significant_umi_records.csv")
        tenx.count_umis_per_cell(
            path,
            meta,
            bam2fasta_args.CELL_BARCODE_PATTERN,
            bam2fasta_args.MOLECULAR_BARCODE_PATTERN,
            3,
            good_barcodes)
        expected_meta = [6, 15, 2, 2, 5, 4, 6, 2]
        all_barcodes = [
            'AAAGATGCAGATCTGT-1', 'AAATGCCCAAACTGCT-1', 'AACACGTAGTGTACCT-1',
            'AACCATGAGTTGTCGT-1', 'AAATGCCGTGAACCTT-1', 'AAATGCCAGATAGTCA-1',
            'AAACGGGAGGATATAC-1', 'AAACGGGTCTCGTATT-1']
        assert all_barcodes == pd.read_csv(
            meta, header=None).iloc[:, 0].values.tolist()
        assert expected_meta == pd.read_csv(
            meta, header=None).iloc[:, 1].values.tolist()
        assert expected_good_barcodes == pd.read_csv(
            good_barcodes, header=None).iloc[:, 0].values.tolist()


def test_record_to_fastq_string():
    path = utils.get_test_data('10x-example/possorted_genome_bam.fastq.gz')
    expected = (
        "@A00111:50:H2H5YDMXX:2:1334:1886:36072\tUB:Z:AGAAAATACG\n" +
        "GATTACTTAGTAGCTGTTTACTTAGCAGCACATTTGCAACAGCATCAAAAGCTATGT" +
        "TACTATAAAATCAGTGCGTGAAGTCTGATTTAC\n")
    with screed.open(path) as f:
        for record_count, record in enumerate(f):
            result = tenx.record_to_fastq_string(record)
            assert expected in result
            break


def test_write_fastq():
    path = utils.get_test_data('10x-example/possorted_genome_bam.fastq.gz')
    with utils.TempDirectory() as location:
        with screed.open(path) as f:
            records = []
            for record_count, record in enumerate(f):
                records.append(record)
        write_path = os.path.join(location, "result.fastq")
        tenx.write_fastq(records, write_path)
        with screed.open(write_path) as f:
            records_written = []
            for record_count, record in enumerate(f):
                records_written.append(record)
            assert records_written == records


def test_make_per_cell_fastqs():
    path = utils.get_test_data('10x-example/possorted_genome_bam.fastq.gz')
    barcodes_file = utils.get_test_data('10x-example/barcodes.tsv')

    barcodes = tenx.read_barcodes_file(barcodes_file)
    with utils.TempDirectory() as location:
        outdir = os.path.join(location, "outdir")
        os.makedirs(outdir)
        tenx.make_per_cell_fastqs(
            path,
            outdir,
            "",
            "fastq",
            bam2fasta_args.CELL_BARCODE_PATTERN,
            barcodes_file)
        fastas = glob.glob(os.path.join(outdir, "*.fastq"))

        for fasta in fastas:
            assert os.path.basename(
                fasta).replace(".fastq", "").replace("_", "") in barcodes


def test_make_per_cell_fastq_gzs():
    path = utils.get_test_data('10x-example/possorted_genome_bam.fastq.gz')
    barcodes_file = utils.get_test_data('10x-example/barcodes.tsv')

    barcodes = tenx.read_barcodes_file(barcodes_file)
    with utils.TempDirectory() as location:
        outdir = os.path.join(location, "outdir")
        os.makedirs(outdir)
        tenx.make_per_cell_fastqs(
            path,
            outdir,
            "possorted_aligned_",
            "fastq.gz",
            bam2fasta_args.CELL_BARCODE_PATTERN,
            barcodes_file)
        fastas = glob.glob(os.path.join(outdir, "*.fastq.gz"))
        for fasta in fastas:
            fasta_name = os.path.basename(fasta).replace(
                ".fastq.gz", "").replace("possorted_aligned__", "")
            assert fasta_name in barcodes
