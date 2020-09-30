import glob
import os
import shutil

import screed

from bam2fasta.tests import bam2fasta_tst_utils as utils
from bam2fasta import cli
from bam2fasta import VERSION


def test_bam2fasta_info():
    status, out, err = utils.run_shell_cmd('bam2fasta info')

    assert "bam2fasta version" in out
    assert "loaded from path" in out
    assert VERSION in out


def test_bam2fasta_info_verbose():
    status, out, err = utils.run_shell_cmd('bam2fasta info -v')

    assert "pathos version" in out
    assert "screed version" in out
    assert "pysam version" in out
    assert "loaded from path" in out


def test_run_count_umis_percell():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data(
            '10x-example/possorted_genome_bam.fastq.gz')
        csv_path = os.path.join(location, "all_barcodes_meta.csv")
        good_barcodes_path = os.path.join(location, "good_barcodes.csv")

        barcodes_path = utils.get_test_data('10x-example/barcodes.tsv')
        renamer_path = utils.get_test_data('10x-example/barcodes_renamer.tsv')

        status, out, err = utils.run_shell_cmd(
            'bam2fasta count_umis_percell --filename ' + testdata1 +
            ' --min-umi-per-barcode 10' +
            ' --write-barcode-meta-csv ' + csv_path +
            ' --barcodes-file ' + barcodes_path + ' --rename-10x-barcodes ' +
            renamer_path + " --processes 1 " +
            "--barcodes-significant-umis-file " + good_barcodes_path,
            in_directory=location)
        assert status == 0
        with open(csv_path, 'rb') as f:
            data = [line.split() for line in f]
        assert len(data) == 8
        with open(good_barcodes_path, 'rb') as f:
            data = [line.split() for line in f]
        assert len(data) == 1


def test_run_make_fastqs_percell():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data(
            '10x-example/possorted_genome_bam.fastq.gz')
        good_barcodes_path = utils.get_test_data(
            '10x-example/good_barcodes.csv')

        barcodes_path = utils.get_test_data('10x-example/barcodes.tsv')
        fastas_dir = os.path.join(location, "fastas")
        if not os.path.exists(fastas_dir):
            os.makedirs(fastas_dir)

        status, out, err = utils.run_shell_cmd(
            'bam2fasta make_fastqs_percell --filename ' + testdata1 +
            ' --barcodes-file ' + barcodes_path +
            " --barcodes-significant-umis-file " + good_barcodes_path +
            ' --save-fastas ' + fastas_dir,
            in_directory=location)
        assert status == 0
        fastqs = glob.glob(os.path.join(fastas_dir + "/*.fastq"))
        assert len(fastqs) == 1, "fastas_dir is {}".format(fastas_dir)


def test_run_bam2fasta_supply_all_args():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example/possorted_genome_bam.bam')
        csv_path = os.path.join(location, "all_barcodes_meta.csv")
        barcodes_path = utils.get_test_data('10x-example/barcodes.tsv')
        renamer_path = utils.get_test_data('10x-example/barcodes_renamer.tsv')
        fastas_dir = os.path.join(location, "fastas")
        temp_fastas_dir = os.path.join(
            os.path.dirname(testdata1), "temp_fastas/")
        if not os.path.exists(fastas_dir):
            os.makedirs(fastas_dir)
        if not os.path.exists(temp_fastas_dir):
            os.makedirs(temp_fastas_dir)

        status, out, err = utils.run_shell_cmd(
            'bam2fasta percell --filename ' + testdata1 +
            ' --min-umi-per-barcode 10' +
            ' --write-barcode-meta-csv ' + csv_path +
            ' --save-intermediate-files ' + temp_fastas_dir +
            ' --barcodes-file ' + barcodes_path + ' --rename-10x-barcodes ' +
            renamer_path + ' --save-fastas ' + fastas_dir + " --processes 1",
            in_directory=location)

        assert status == 0
        with open(csv_path, 'rb') as f:
            data = [line.split() for line in f]
        assert len(data) == 9
        fasta_files = os.listdir(fastas_dir)
        barcodes = [filename.replace(".fasta", "") for filename in fasta_files]
        assert len(barcodes) == 1
        assert len(fasta_files) == 1
        assert barcodes[0] == \
            ('lung_epithelial_cell_AAATGCCCAAACTGCT-1_bam2fasta')
        count = 0
        fasta_file_name = os.path.join(fastas_dir, fasta_files[0])
        for record in screed.open(fasta_file_name):
            name = record.name
            sequence = record.sequence
            count += 1
            assert name.startswith('lung_epithelial_cell_AAATGCCCAAACTGCT-1')
            assert sequence.count(">") == 0
            assert sequence.count("X") == 0
        shutil.rmtree(temp_fastas_dir)


def test_run_bam2fasta_default_args():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example/possorted_genome_bam.bam')

        status, out, err = utils.run_shell_cmd(
            'bam2fasta percell --filename ' + testdata1,
            in_directory=location)

        assert status == 0
        fasta_files = os.listdir(location)
        barcodes = [
            filename.replace(".fasta", "") for
            filename in fasta_files if filename.endswith("_bam2fasta.fasta")]
        assert len(barcodes) == 8


def test_run_bam2fasta_percell():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example/possorted_genome_bam.bam')

        fasta_files = cli.percell(
            ['--filename', testdata1, '--save-fastas', location,
             ])

        barcodes = [filename.replace(".fasta", "") for filename in fasta_files]
        assert len(barcodes) == 8


def test_run_bam2fasta_convert():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example/possorted_genome_bam.bam')

        fasta_files = cli.convert(
            ['--filename', testdata1, '--save-fastas', location,
             ])

        barcodes = [filename.replace(".fasta", "") for filename in fasta_files]
        assert len(barcodes) == 8


def test_run_bam2fasta_percell_no_shard():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data(
            '10x-example/possorted_genome_bam.fastq.gz')

        fasta_files = cli.percell(
            ['--filename', testdata1, '--save-fastas', location])
        print(fasta_files)

        barcodes = [filename.replace(".fasta", "") for filename in fasta_files]
        assert len(barcodes) == 8


def test_run_bam2fasta_percell_nonzero_umi():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example/possorted_genome_bam.bam')

        fasta_files = cli.percell(
            ['--filename', testdata1, '--save-fastas', location,
             '--min-umi-per-barcode', '10', ])
        barcodes = [filename.replace(".fasta", "") for filename in fasta_files]
        assert len(barcodes) == 1
        sequences_fasta = []
        with screed.open(fasta_files[0]) as f:
            for record in f:
                sequences_fasta.append(record.sequence)
        gt_data = utils.get_test_data(
            '10x-example/groundtruth_fasta_sequences.txt')
        with open(gt_data, "r") as f:
            for index, line in enumerate(f.readlines()):
                assert line.strip() in sequences_fasta, \
                    "failed at index {}".format(index)


def test_run_bam2fasta_fq_percell_no_shard_nonzero_umi():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data(
            '10x-example/possorted_genome_bam.fastq.gz')

        fasta_files = cli.percell(
            ['--filename', testdata1, '--save-fastas', location,
             '--min-umi-per-barcode', '10'])
        print(fasta_files)
        barcodes = [
            filename.replace(".fastq", "") for filename in fasta_files]
        assert len(barcodes) == 1
        sequences_fastq = []
        with screed.open(fasta_files[0]) as f:
            for record in f:
                sequences_fastq.append(record.sequence)
        gt_data = utils.get_test_data(
            '10x-example/groundtruth_fasta_sequences.txt')
        with open(gt_data, "r") as f:
            for index, line in enumerate(f.readlines()):
                assert line.strip() in sequences_fastq, \
                    "failed at index {}".format(index)


def test_run_bam2fasta_fq_gz_percell_no_shard_nonzero_umi():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data(
            '10x-example/possorted_genome_bam.fastq.gz')

        fasta_files = cli.percell(
            ['--filename', testdata1, '--save-fastas', location,
             '--min-umi-per-barcode', '10', '--output-format', 'fastq.gz'])
        barcodes = [
            filename.replace(".fastq.gz", "") for filename in fasta_files]
        assert len(barcodes) == 1
        sequences_fastq = []
        with screed.open(fasta_files[0]) as f:
            for record in f:
                sequences_fastq.append(record.sequence)
        gt_data = utils.get_test_data(
            '10x-example/groundtruth_fasta_sequences.txt')
        with open(gt_data, "r") as f:
            for index, line in enumerate(f.readlines()):
                assert line.strip() in sequences_fastq, \
                    "failed at index {}".format(index)
